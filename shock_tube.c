/*
 * shock_tube.c — 1D Shock Tube solver
 *
 * Two modes:
 *   ./shock_tube        — MUSCL-HLLC with SSP-RK2  (default)
 *   ./shock_tube --lf   — Lax-Friedrich with forward Euler
 *
 * Both use the same real-gas thermodynamic model (tables + compression factor).
 * Output: profiles[_lf].dat, timeseries[_lf].dat
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "tables.h"

/* ========== Physical constants ========== */
#define CV    1647.0
#define CP    2160.0
#define RGAS  (8314.46 / 16.38)
#define GAMMA (CP / CV)

/* Ideal-gas consistent gamma: enforces p = (GAMMA_I-1)*rho*e_int = rho*RGAS*T
   with e_int = CV*T.  For a real gas CP-CV != R, so GAMMA != 1+R/CV. */
#define GAMMA_I (1.0 + RGAS / CV)

#define PI_VAL    3.14159265358979323846

/* Wall thermal model: 1D radial conduction through steel + concrete */
#define N_WALL_STEEL 5
#define N_WALL_CONC  3
#define N_WALL       (N_WALL_STEEL + N_WALL_CONC)
#define K_STEEL      50.0     /* W/(m·K) — thermal conductivity */
#define RHO_STEEL    7800.0   /* kg/m³ */
#define CP_STEEL     500.0    /* J/(kg·K) */
#define T_STEEL_WALL 0.0268   /* m (26.8 mm wall thickness at rupture sites) */
#define K_CONC       1.5      /* W/(m·K) — concrete weight coating */
#define RHO_CONC     2400.0   /* kg/m³ */
#define CP_CONC      880.0    /* J/(kg·K) */
#define T_CONC_WALL  0.100    /* m (100 mm reinforced concrete coating) */

/* Gas thermal properties for convective HTC (Dittus-Boelter) */
#define K_GAS_CONV   0.045    /* W/(m·K) — approx for CH4 at moderate pressure */
#define MU_GAS       1.3e-5   /* Pa·s — dynamic viscosity */
#define PR_GAS       0.8      /* Prandtl number */

/* Numerics */
#define CFL_NUM    0.6

/* Output */
#define N_TSERIES  5000
static int g_n_profiles = 100;

/* Maximum number of cells */
#define MAX_CELLS 500000

/* ========== Runtime parameters (set from CLI or defaults) ========== */

/* Pipe geometry & conditions (parameterized via CLI) */
static double g_dia         = 1.153;
static double g_rough       = 0.017e-3;
static double g_pressure_L  = 165e5;
static double g_pressure_R  = 7.1e5;
static double g_T0          = 282.0;

/* Derived geometry (computed once after CLI parsing) */
static double g_area;
static double g_circum;

/* Wall temperature array: N_WALL radial nodes per gas cell, row-major */
static double g_Tw[MAX_CELLS * N_WALL];
#define TW(i,j) g_Tw[(i)*N_WALL + (j)]

/* Grid parameters */
static double g_pipe_length  = 4400.0;
static double g_end_time     = 20.0;
static double g_dx_init      = 1.0;
static double g_dx_init_left = 0.0;  /* 0 = one-sided grid (default) */
static double g_dx_max       = 0.5;
static double g_fac          = 1.0;
static double g_dx_trans     = 1000.0;  /* transition length from dx_init to dx_max [m] */

/* Global flag to disable source terms for diagnostics */
static int g_no_source = 0;
/* Global flag to disable heat transfer only (keep friction) */
static int g_no_heat = 0;
/* Global flag: use old simplified constant-T_wall heat transfer model */
static int g_simple_wall = 0;
/* Global flag for ideal gas mode (bypass tables, Z=1) */
static int g_ideal_gas = 0;
/* Global flag to disable choked outlet BC */
static int g_no_choke = 0;
/* Override friction factor from CLI (0 = use computed value) */
static double g_lambda_override = 0.0;
/* Diagnostic: track choked BC state */
static double g_poutlet_choked = 0.0;
static double g_sim_time = 0.0;  /* accessible from BC function for diagnostics */
static int g_choked_active = 0;

/* Adaptive merge parameters */
static int    g_merge        = 0;      /* enable adaptive merging */
static double g_merge_tol    = 0.01;   /* relative tolerance for merging */
static double g_merge_dx_max = 5000.0; /* max cell size after merge (m) */
static int    g_i_front      = 0;      /* leftmost active cell index */

/* Early output parameters */
static double g_early_end = 0.0;  /* 0 = disabled */
static double g_early_dt  = 0.01; /* profile interval during early phase */

/* Stopping criterion: stop when closed-end pressure drops below this (Pa). 0=disabled */
static double g_pstop = 0.0;

/* Tube mode: extend pipe with vacuum section (classic shock tube problem).
   When > 0, domain goes from -pipe_length to +vacuum_ext.
   Both ends are reflective walls, diaphragm at x=0 opens at t=0. */
static double g_vacuum_ext = 0.0;
static int    g_i_diaphragm = -1;  /* cell index nearest x=0 */
static double g_tube_start  = 0.0; /* >0: tube mode for this many seconds, then switch to outlet BC */

/* CoolProp internal energy offset: U_coolprop(T) ≈ CV*T + g_U_offset at low P.
   Computed at startup for ideal-gas fallback when T is outside table range. */
static double g_U_offset = 0.0;

/* ========== Table lookup: bicubic Catmull-Rom interpolation (C1 smooth) ========== */

/* Catmull-Rom basis: given 4 values f[-1], f[0], f[1], f[2] and fraction t in [0,1],
   returns interpolated value at f[0] + t*(f[1]-f[0]) with smooth derivatives. */
static inline double catmull_rom(double fm1, double f0, double f1, double f2, double t)
{
    double t2 = t * t, t3 = t2 * t;
    return 0.5 * ((2.0 * f0) +
                   (-fm1 + f1) * t +
                   (2.0 * fm1 - 5.0 * f0 + 4.0 * f1 - f2) * t2 +
                   (-fm1 + 3.0 * f0 - 3.0 * f1 + f2) * t3);
}

static inline double table_lookup(const double tab[NT][NP], double T, double P)
{
    /* Clamp to grid bounds */
    if (T < T_MIN) T = T_MIN;
    if (T > T_MIN + (NT - 1) * T_STEP) T = T_MIN + (NT - 1) * T_STEP;
    if (P < P_MIN) P = P_MIN;
    if (P > P_MIN + (NP - 1) * P_STEP) P = P_MIN + (NP - 1) * P_STEP;

    double fi = (T - T_MIN) / T_STEP;
    double fj = (P - P_MIN) / P_STEP;

    int i = (int)fi;
    int j = (int)fj;

    if (i >= NT - 1) i = NT - 2;
    if (j >= NP - 1) j = NP - 2;

    double wt = fi - i;
    double wp = fj - j;

    /* Clamp indices for the 4x4 stencil to stay within table bounds */
    int im1 = (i > 0) ? i - 1 : 0;
    int ip2 = (i < NT - 2) ? i + 2 : NT - 1;
    int jm1 = (j > 0) ? j - 1 : 0;
    int jp2 = (j < NP - 2) ? j + 2 : NP - 1;

    /* Interpolate 4 rows in P direction, then interpolate results in T direction */
    double r0 = catmull_rom(tab[im1][jm1], tab[im1][j], tab[im1][j+1], tab[im1][jp2], wp);
    double r1 = catmull_rom(tab[i  ][jm1], tab[i  ][j], tab[i  ][j+1], tab[i  ][jp2], wp);
    double r2 = catmull_rom(tab[i+1][jm1], tab[i+1][j], tab[i+1][j+1], tab[i+1][jp2], wp);
    double r3 = catmull_rom(tab[ip2][jm1], tab[ip2][j], tab[ip2][j+1], tab[ip2][jp2], wp);

    return catmull_rom(r0, r1, r2, r3, wt);
}

static inline double lookup_U(double T, double P)  {
    if (g_ideal_gas) return CV * T;
    return table_lookup(U_table, T, P);
}
static inline double lookup_Cv(double T, double P)  {
    (void)P;
    if (g_ideal_gas) return CV;
    double cv = table_lookup(Cv_table, T, P);
    if (cv < 500.0) cv = 500.0;  /* safety floor near critical region */
    return cv;
}
/* ========== Compression factor (analytical, 19 coefficients) ========== */

static inline double compression_factor(double T, double P)
{
    if (g_ideal_gas) { (void)T; (void)P; return 1.0; }
    if (T < 150.0) T = 150.0;
    if (P < 100.0) P = 100.0;

    const double gammag = 0.5654;
    const double Pc = 756.8 - 131.0 * gammag - 3.6 * gammag * gammag;
    const double Tc = 169.2 + 349.5 * gammag - 74.0 * gammag * gammag;

    const double a1  =  0.317842;
    const double a2  =  0.382216;
    const double a3  = -7.76835;
    const double a4  = 14.2905;
    const double a5  =  2e-6;
    const double a6  = -0.004693;
    const double a7  =  0.096254;
    const double a8  =  0.16672;
    const double a9  =  0.96691;
    const double a10 =  0.063069;
    const double a11 = -1.96685;
    const double a12 = 21.0581;
    const double a13 = -27.0246;
    const double a14 = 16.23;
    const double a15 = 207.783;
    const double a16 = -488.161;
    const double a17 = 176.29;
    const double a18 =  1.88453;
    const double a19 =  3.05921;

    double Tpr = (T * 1.8) / Tc;
    double Ppr = (P / 100000.0 * 14.5) / Pc;

    double t = 1.0 / Tpr;

    double A = a1 * exp(a2 * (1.0 - t) * (1.0 - t)) * Ppr * t;
    double B = a3 * t + a4 * t * t + a5 * pow(Ppr, 6.0) * pow(t, 6.0);
    double Cc = a9 + a8 * Ppr * t + a7 * Ppr * Ppr * t * t + a6 * Ppr * Ppr * Ppr * t * t * t;
    double Dc = a10 * exp(a11 * (1.0 - t) * (1.0 - t)) * t;
    double Ec = a12 * t + a13 * t * t + a14 * t * t * t;
    double F  = a15 * t + a16 * t * t + a17 * t * t * t;
    double G  = a18 + a19 * t;

    double y = (Dc * Ppr) / (-(A * A * B) / (Cc * Cc * Cc) + (1.0 + A * A) / Cc);

    double Z = (Dc * Ppr * (1.0 + y + y * y - y * y * y))
             / (pow(1.0 - y, 3.0) * (Dc * Ppr + Ec * y * y - F * pow(y, G)));

    return Z;
}

/* ========== Thermodynamic sound speed from EOS ========== */
/*
 * For p = rho * Z(T,P) * R * T:
 *   D = 1 - p * (dZ/dP)_T / Z
 *   (dp/drho)_T = Z * R * T / D
 *   (dp/dT)_rho = rho * R * (T * (dZ/dT)_P + Z) / D
 *   a^2 = (dp/drho)_T + T * [(dp/dT)_rho]^2 / (rho^2 * cv)
 */
static inline double thermo_soundspeed(double T, double P, double rho_val, double cv)
{
    if (g_ideal_gas) return sqrt(GAMMA_I * P / rho_val);

    /* At very low pressure, use ideal gas (Z≈1, accurate) */
    if (P < 1000.0) return sqrt(GAMMA_I * P / rho_val);

    double Z = compression_factor(T, P);

    /* Central differences for dZ/dT and dZ/dP */
    double hT = 0.02;   /* K */
    double hP = 200.0;  /* Pa */
    double dZdT = (compression_factor(T + hT, P) - compression_factor(T - hT, P)) / (2.0 * hT);
    double dZdP = (compression_factor(T, P + hP) - compression_factor(T, P - hP)) / (2.0 * hP);

    double D = 1.0 - P * dZdP / Z;
    if (fabs(D) < 1e-10) D = 1e-10;

    double dpdrho_T = Z * RGAS * T / D;
    double dpdT_rho = rho_val * RGAS * (T * dZdT + Z) / D;

    double a2 = dpdrho_T + T * dpdT_rho * dpdT_rho / (rho_val * rho_val * cv);

    if (a2 < 0.0) a2 = GAMMA * P / rho_val;  /* fallback */
    return sqrt(a2);
}

/* ========== Van Leer limiter (smoother than MC) ========== */

static inline double vanleer_limiter(double a, double b)
{
    if (a * b <= 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
}

/* ========== Darcy friction factor ========== */

static double compute_lambda(void)
{
    double x = -2.0 * log10(g_rough / (3.7 * g_dia));
    return 1.0 / (x * x);
}

/* ========== Global arrays ========== */

static double x_face[MAX_CELLS + 1];  /* face positions */
static double xc[MAX_CELLS];          /* cell centres */
static double dx[MAX_CELLS];          /* cell widths */

static double rho[MAX_CELLS], u[MAX_CELLS], e[MAX_CELLS];
static double T_arr[MAX_CELLS], p[MAX_CELLS], a_arr[MAX_CELLS];

/* Work arrays for RK stages */
static double rho1[MAX_CELLS], u1[MAX_CELLS], e1[MAX_CELLS];
static double T1[MAX_CELLS], p1_arr[MAX_CELLS], a1_arr[MAX_CELLS];

static double rho2[MAX_CELLS], u2_arr[MAX_CELLS], e2[MAX_CELLS];

static double rho_0[MAX_CELLS], mom_0[MAX_CELLS], e_0[MAX_CELLS];

/* Flux arrays */
static double rho_flux[MAX_CELLS], mom_flux[MAX_CELLS], e_flux[MAX_CELLS];

/* Slope arrays (conserved reconstruction: rho, mom, e) */
static double s_rho[MAX_CELLS], s_mom[MAX_CELLS], s_e[MAX_CELLS];

/* Cv saved for heat transfer */
static double Cv_real[MAX_CELLS], Cv_real1[MAX_CELLS];

/* ========== Outlet BC: choked/subsonic isentropic exit ========== */

/* Computes the outlet face flux using isentropic relations.
   When the sonic pressure exceeds back-pressure → choked (M=1) exit.
   When sonic pressure < back-pressure → subsonic exit at back-pressure.
   This gives a smooth, continuous transition at the choking boundary. */
/* Set ghost cell N-1 to the proper outlet BC state.
   When choked: ghost cell is the sonic isentropic exit state (M=1).
   When subsonic: ghost cell is at back-pressure p_R.
   The HLLC/LF naturally computes the correct outlet face flux from this. */
static void set_outlet_ghost(double *rho_arr, double *u_arr, double *e_arr,
                             double *p_arr, double *T_arr, double *a_arr,
                             double *Cv_arr, int N)
{
    /* Tube mode: reflective wall at right end (mirror of N-2, negate velocity) */
    if (g_vacuum_ext > 0.0) {
        rho_arr[N-1] = rho_arr[N-2];
        u_arr[N-1]   = -u_arr[N-2];
        p_arr[N-1]   = p_arr[N-2];
        T_arr[N-1]   = T_arr[N-2];
        a_arr[N-1]   = a_arr[N-2];
        Cv_arr[N-1]  = Cv_arr[N-2];
        e_arr[N-1]   = e_arr[N-2];  /* same total energy (KE same due to |u|) */
        g_choked_active = 0;
        g_poutlet_choked = 0.0;
        return;
    }

    double gam = g_ideal_gas ? GAMMA_I : GAMMA;
    double gm1 = gam - 1.0;

    double u_approach = fmax(u_arr[N-2], 0.0);
    /* Use thermodynamic sound speed (consistent with solver) */
    double a_cell = a_arr[N-2];
    double M_cell = u_approach / a_cell;
    if (M_cell > 0.999) M_cell = 0.999;

    double M2 = M_cell * M_cell;
    double stag_fac = 1.0 + 0.5 * gm1 * M2;
    double p0 = p_arr[N-2] * pow(stag_fac, gam / gm1);
    double p_sonic = p0 * pow(2.0 / (gam + 1.0), gam / gm1);
    g_poutlet_choked = p_sonic;

    if (!g_no_choke && p_sonic > g_pressure_R) {
        /* Choked: ghost cell at isentropic sonic exit state (M=1). */
        g_choked_active = 1;
        double T_val = g_ideal_gas ? p_arr[N-2] / (rho_arr[N-2] * RGAS) : T_arr[N-2];
        double T_stag = T_val * stag_fac;
        double T_exit = T_stag * 2.0 / (gam + 1.0);
        double p_exit = p_sonic;
        double a_exit;

        if (g_ideal_gas) {
            a_exit = sqrt(gam * RGAS * T_exit);
            rho_arr[N-1] = p_exit / (RGAS * T_exit);
            Cv_arr[N-1] = CV;
        } else {
            double T_c = fmax(T_exit, 150.0);
            double P_c = fmax(p_exit, 1e5);
            double cfac = compression_factor(T_c, P_c);
            rho_arr[N-1] = p_exit / (cfac * RGAS * T_c);
            Cv_arr[N-1] = lookup_Cv(T_c, P_c);
            a_exit = thermo_soundspeed(T_c, P_c, rho_arr[N-1], Cv_arr[N-1]);
        }

        u_arr[N-1] = a_exit;
        p_arr[N-1] = p_exit;
        T_arr[N-1] = T_exit;
        a_arr[N-1] = a_exit;
        e_arr[N-1] = rho_arr[N-1] * (g_ideal_gas ? CV * T_exit : lookup_U(fmax(T_exit,150.0), fmax(p_exit,1e5)))
                   + 0.5 * rho_arr[N-1] * a_exit * a_exit;
    } else {
        /* Subsonic: characteristic-based (NRBC) exit at back-pressure p_R.
           Extrapolate outgoing Riemann invariant J+ = u + 2a/(gm1) from interior,
           compute boundary sound speed at p_R using isentropic relation,
           then u_exit = J+ - 2*a_exit/(gm1).
           This is smooth at the choked→subsonic transition because at M=1,
           J+ = a*(1 + 2/(gm1)) = a*(gam+1)/(gm1), giving u = a. */
        g_choked_active = 0;

        /* Outgoing Riemann invariant from cell N-2 */
        double Jplus = fmax(u_arr[N-2], 0.0) + 2.0 * a_arr[N-2] / gm1;

        /* Isentropic sound speed at p_R: a_exit/a_int = (p_R/p_int)^((gm1)/(2*gam)) */
        double p_ratio = g_pressure_R / fmax(p_arr[N-2], 1e5);
        double a_exit = a_arr[N-2] * pow(p_ratio, gm1 / (2.0 * gam));

        /* Exit velocity from Riemann invariant */
        double u_exit = Jplus - 2.0 * a_exit / gm1;

        /* Isentropic temperature at p_R */
        double T_val = g_ideal_gas ? p_arr[N-2] / (rho_arr[N-2] * RGAS) : T_arr[N-2];
        double T_c = T_val * pow(p_ratio, gm1 / gam);
        if (T_c < 150.0) T_c = 150.0;
        if (T_c > 400.0) T_c = 400.0;

        if (g_ideal_gas) {
            a_exit = sqrt(gam * RGAS * T_c);
            rho_arr[N-1] = g_pressure_R / (RGAS * T_c);
            Cv_arr[N-1] = CV;
        } else {
            double P_c = fmax(g_pressure_R, 1e5);
            double cfac_bc = compression_factor(T_c, P_c);
            rho_arr[N-1] = g_pressure_R / (cfac_bc * RGAS * T_c);
            Cv_arr[N-1] = lookup_Cv(T_c, P_c);
            a_exit = thermo_soundspeed(T_c, P_c, rho_arr[N-1], Cv_arr[N-1]);
            /* Recompute u_exit with the accurate real-gas sound speed */
            u_exit = Jplus - 2.0 * a_exit / gm1;
        }

        u_arr[N-1] = u_exit;
        p_arr[N-1] = g_pressure_R;
        T_arr[N-1] = T_c;
        a_arr[N-1] = a_exit;
        e_arr[N-1] = rho_arr[N-1] * (g_ideal_gas ? CV * T_c : lookup_U(T_c, fmax(g_pressure_R, 1e5)))
                   + 0.5 * rho_arr[N-1] * u_exit * u_exit;
    }
}

/* ========== Wall conduction + heat exchange ========== */

/*
 * 1D radial conduction through composite wall (steel + concrete).
 * Inner BC: Dirichlet T_surface = T_gas (Bi >> 1, gas-side convection dominates).
 * Outer BC: insulated (seawater convection negligible on relevant timescales).
 *
 * Node layout (per gas cell):
 *   j = 0..N_WALL_STEEL-1 : steel, spacing dr_s = T_STEEL_WALL / N_WALL_STEEL
 *   j = N_WALL_STEEL..N_WALL-1 : concrete, spacing dr_c = T_CONC_WALL / N_WALL_CONC
 *
 * Heat flux to gas = k_steel * (TW[0] - T_gas) / (dr_s/2) per unit inner area.
 */
static void update_wall_and_heat(int N, int i_start, double *e_arr,
                                 const double *T_gas, const double *Cv_arr,
                                 const double *rho_arr, const double *u_arr,
                                 double dt)
{
    if (g_no_source || g_no_heat) return;

    const double circ_over_area = g_circum / g_area;

    /* --- Simplified constant-T_wall model (old model) --- */
    if (g_simple_wall) {
        /* Steady-state conduction through steel+concrete wall to seawater at T0.
         * Total thermal resistance: steel/k_steel + concrete/k_conc [m²·K/W] */
        const double R_wall = T_STEEL_WALL / K_STEEL + T_CONC_WALL / K_CONC;
        const double U_wall = (1.0 / R_wall) * circ_over_area;  /* W/(m³·K) */
        for (int i = i_start; i < N - 1; i++) {
            double tdiff = g_T0 - T_gas[i];
            double tau = U_wall / (rho_arr[i] * Cv_arr[i]);
            double deltaT = tdiff * (1.0 - exp(-dt * tau));
            e_arr[i] += deltaT * rho_arr[i] * Cv_arr[i];
        }
        return;
    }

    /* --- Full transient wall model: 1D radial conduction --- */
    const int ns = N_WALL_STEEL;
    const double dr_s = T_STEEL_WALL / N_WALL_STEEL;
    const double dr_c = T_CONC_WALL / N_WALL_CONC;
    const double alpha_s = K_STEEL / (RHO_STEEL * CP_STEEL);
    const double alpha_c = K_CONC / (RHO_CONC * CP_CONC);
    /* Interface thermal resistance: half-cell in steel + half-cell in concrete */
    const double R_intf = dr_s / (2.0 * K_STEEL) + dr_c / (2.0 * K_CONC);
    /* Dittus-Boelter pre-factor: 0.023 * Pr^0.4 * k_gas / D */
    const double DB_coeff = 0.023 * pow(PR_GAS, 0.4) * K_GAS_CONV / g_dia;
    const double half_dr_over_k = 0.5 * dr_s / K_STEEL;  /* conductive resistance */

    for (int i = i_start; i < N - 1; i++) {
        double Tg = T_gas[i];

        /* --- Flow-dependent convective HTC (Dittus-Boelter) --- */
        double Re = rho_arr[i] * fabs(u_arr[i]) * g_dia / MU_GAS;
        double h_conv = DB_coeff * pow(fmax(Re, 1.0), 0.8);
        if (h_conv < 5.0) h_conv = 5.0;  /* natural convection floor */

        /* Robin BC: effective conductance gas → node-0 center [W/(m²·K)] */
        double U_eff = 1.0 / (1.0 / h_conv + half_dr_over_k);

        /* --- Heat flux from wall to gas (using OLD wall state) --- */
        double q_vol = U_eff * (TW(i,0) - Tg) * circ_over_area;  /* W/m³ */
        e_arr[i] += q_vol * dt;

        /* --- Update wall temperatures (explicit 1D conduction) --- */
        double Tnew[N_WALL];

        /* Steel node 0: Robin BC (convective inner surface)
           Energy balance: rho*cp*dr * dT/dt = U_eff*(Tg-T0) + k*(T1-T0)/dr */
        Tnew[0] = TW(i,0) + dt / (RHO_STEEL * CP_STEEL * dr_s)
                   * (U_eff * (Tg - TW(i,0))
                      + K_STEEL * (TW(i,1) - TW(i,0)) / dr_s);

        /* Interior steel nodes */
        for (int j = 1; j < ns - 1; j++)
            Tnew[j] = TW(i,j) + alpha_s * dt / (dr_s * dr_s)
                       * (TW(i,j-1) - 2.0 * TW(i,j) + TW(i,j+1));

        /* Last steel node (interface with concrete) */
        {
            double flux_left  = K_STEEL * (TW(i,ns-2) - TW(i,ns-1)) / dr_s;
            double flux_right = (TW(i,ns-1) - TW(i,ns)) / R_intf;
            Tnew[ns-1] = TW(i,ns-1) + dt / (RHO_STEEL * CP_STEEL * dr_s)
                          * (flux_left - flux_right);
        }

        /* First concrete node (interface with steel) */
        {
            double flux_left  = (TW(i,ns-1) - TW(i,ns)) / R_intf;
            double flux_right = (ns + 1 < N_WALL)
                ? K_CONC * (TW(i,ns) - TW(i,ns+1)) / dr_c : 0.0;
            Tnew[ns] = TW(i,ns) + dt / (RHO_CONC * CP_CONC * dr_c)
                        * (flux_left - flux_right);
        }

        /* Interior concrete nodes */
        for (int j = ns + 1; j < N_WALL - 1; j++)
            Tnew[j] = TW(i,j) + alpha_c * dt / (dr_c * dr_c)
                       * (TW(i,j-1) - 2.0 * TW(i,j) + TW(i,j+1));

        /* Last concrete node: Dirichlet outer surface (seawater at T0) */
        Tnew[N_WALL-1] = TW(i,N_WALL-1) + alpha_c * dt / (dr_c * dr_c)
                          * (TW(i,N_WALL-2) - 2.0 * TW(i,N_WALL-1) + g_T0);

        /* Copy back */
        for (int j = 0; j < N_WALL; j++) TW(i,j) = Tnew[j];
    }
}

/* ========== HLLC flux computation ========== */

static void compute_hllc_flux(
    const double *rho_in, const double *u_in, const double *e_in,
    const double *p_in,
    const double *T_in, const double *Cv_in,
    int N, int i_start, double dt, double lambda_fric,
    int is_stage2 __attribute__((unused)),
    /* outputs: */
    double *rho_out, double *u_out, double *e_out)
{
    double mom[MAX_CELLS];
    for (int i = i_start; i < N; i++) mom[i] = rho_in[i] * u_in[i];

    /* MUSCL reconstruction with van Leer limiter on CONSERVED variables */
    s_rho[i_start] = 0.0; s_mom[i_start] = 0.0; s_e[i_start] = 0.0;
    s_rho[N-1] = 0.0; s_mom[N-1] = 0.0; s_e[N-1] = 0.0;
    for (int i = i_start + 1; i < N - 1; i++) {
        double dxL = xc[i] - xc[i-1];
        double dxR = xc[i+1] - xc[i];
        s_rho[i] = vanleer_limiter((rho_in[i] - rho_in[i-1]) / dxL, (rho_in[i+1] - rho_in[i]) / dxR);
        s_mom[i] = vanleer_limiter((mom[i] - mom[i-1]) / dxL, (mom[i+1] - mom[i]) / dxR);
        s_e[i]   = vanleer_limiter((e_in[i] - e_in[i-1]) / dxL, (e_in[i+1] - e_in[i]) / dxR);
    }

    /* Compute HLLC fluxes at each face */
    int Nfaces = N - 1;
    for (int f = i_start; f < Nfaces; f++) {
        int iL = f, iR = f + 1;

        /* Reconstruct conserved left/right states at face */
        double rhoL = rho_in[iL] + 0.5 * dx[iL] * s_rho[iL];
        double rhoR = rho_in[iR] - 0.5 * dx[iR] * s_rho[iR];
        double momL = mom[iL]    + 0.5 * dx[iL] * s_mom[iL];
        double momR = mom[iR]    - 0.5 * dx[iR] * s_mom[iR];
        double eL   = e_in[iL]   + 0.5 * dx[iL] * s_e[iL];
        double eR   = e_in[iR]   - 0.5 * dx[iR] * s_e[iR];

        /* Positivity */
        if (rhoL < 1e-6) rhoL = 1e-6;
        if (rhoR < 1e-6) rhoR = 1e-6;
        if (eL < 1e-6) eL = 1e-6;
        if (eR < 1e-6) eR = 1e-6;

        double uL = momL / rhoL;
        double uR = momR / rhoR;

        /* Recover face pressure and sound speed from reconstructed conserved state */
        double e_int_L = eL / rhoL - 0.5 * uL * uL;
        double e_int_R = eR / rhoR - 0.5 * uR * uR;
        if (e_int_L < 1000.0) e_int_L = 1000.0;
        if (e_int_R < 1000.0) e_int_R = 1000.0;

        double pL, pR, aL, aR;

        if (g_ideal_gas) {
            /* Direct: p = rho*R*T = rho*R*(e_int/Cv) = (GAMMA_I-1)*rho*e_int */
            pL = (GAMMA_I - 1.0) * rhoL * e_int_L;
            pR = (GAMMA_I - 1.0) * rhoR * e_int_R;
            aL = sqrt(GAMMA_I * pL / rhoL);
            aR = sqrt(GAMMA_I * pR / rhoR);
        } else {
            /* Iterated face recovery: Newton for T + update P each iteration.
               Falls back to ideal gas if T wanders outside table range. */
            double TfL = T_in[iL];
            double PfL = fmax(p_in[iL], 1e5);
            for (int it = 0; it < 3; it++) {
                TfL += (e_int_L - lookup_U(TfL, PfL)) / lookup_Cv(TfL, PfL);
                if (TfL < 150.0) TfL = 150.0;
                if (TfL > 400.0) TfL = 400.0;
                PfL = rhoL * compression_factor(TfL, PfL) * RGAS * TfL;
            }
            if (TfL <= 155.0 || TfL >= 395.0) {
                /* Outside table range — use ideal gas with U offset */
                double e_int_IG = e_int_L - g_U_offset;
                if (e_int_IG < 1000.0) e_int_IG = 1000.0;
                pL = (GAMMA_I - 1.0) * rhoL * e_int_IG;
                aL = sqrt(GAMMA_I * pL / rhoL);
            } else {
                pL = PfL;
                aL = thermo_soundspeed(TfL, PfL, rhoL, lookup_Cv(TfL, PfL));
            }

            double TfR = T_in[iR];
            double PfR = fmax(p_in[iR], 1e5);
            for (int it = 0; it < 3; it++) {
                TfR += (e_int_R - lookup_U(TfR, PfR)) / lookup_Cv(TfR, PfR);
                if (TfR < 150.0) TfR = 150.0;
                if (TfR > 400.0) TfR = 400.0;
                PfR = rhoR * compression_factor(TfR, PfR) * RGAS * TfR;
            }
            if (TfR <= 155.0 || TfR >= 395.0) {
                /* Outside table range — use ideal gas with U offset */
                double e_int_IG = e_int_R - g_U_offset;
                if (e_int_IG < 1000.0) e_int_IG = 1000.0;
                pR = (GAMMA_I - 1.0) * rhoR * e_int_IG;
                aR = sqrt(GAMMA_I * pR / rhoR);
            } else {
                pR = PfR;
                aR = thermo_soundspeed(TfR, PfR, rhoR, lookup_Cv(TfR, PfR));
            }
        }

        /* Wave speed estimates */
        double SL = fmin(uL - aL, uR - aR);
        double SR = fmax(uL + aL, uR + aR);

        double denom = rhoL * (SL - uL) - rhoR * (SR - uR);
        if (fabs(denom) < 1e-30) denom = 1e-30;
        double SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / denom;

        /* Physical fluxes */
        double FL_r = momL;
        double FL_m = rhoL * uL * uL + pL;
        double FL_e = uL * (eL + pL);
        double FR_r = momR;
        double FR_m = rhoR * uR * uR + pR;
        double FR_e = uR * (eR + pR);

        double flux_r, flux_m, flux_e;

        if (SL >= 0.0) {
            flux_r = FL_r;
            flux_m = FL_m;
            flux_e = FL_e;
        } else if (SM >= 0.0) {
            /* Star-left */
            double dSL = SL - uL;
            if (fabs(dSL) < 1e-30) dSL = 1e-30;
            double rho_sL = rhoL * (SL - uL) / (SL - SM);
            double mom_sL = rho_sL * SM;
            double E_sL   = rho_sL * (eL / rhoL + (SM - uL) * (SM + pL / (rhoL * dSL)));
            flux_r = FL_r + SL * (rho_sL - rhoL);
            flux_m = FL_m + SL * (mom_sL - momL);
            flux_e = FL_e + SL * (E_sL - eL);
        } else if (SR >= 0.0) {
            /* Star-right */
            double dSR = SR - uR;
            if (fabs(dSR) < 1e-30) dSR = 1e-30;
            double rho_sR = rhoR * (SR - uR) / (SR - SM);
            double mom_sR = rho_sR * SM;
            double E_sR   = rho_sR * (eR / rhoR + (SM - uR) * (SM + pR / (rhoR * dSR)));
            flux_r = FR_r + SR * (rho_sR - rhoR);
            flux_m = FR_m + SR * (mom_sR - momR);
            flux_e = FR_e + SR * (E_sR - eR);
        } else {
            flux_r = FR_r;
            flux_m = FR_m;
            flux_e = FR_e;
        }

        rho_flux[f] = flux_r;
        mom_flux[f] = flux_m;
        e_flux[f]   = flux_e;
    }

    /* Conservative update */
    /* Cell i_start (left boundary: wall/reflective — pressure contributes) */
    rho_out[i_start] = rho_in[i_start] - dt / dx[i_start] * rho_flux[i_start];
    double mom1_val = rho_in[i_start] * u_in[i_start]
                    - dt / dx[i_start] * (mom_flux[i_start] - p_in[i_start]);
    u_out[i_start] = mom1_val / rho_out[i_start];
    e_out[i_start] = e_in[i_start] - dt / dx[i_start] * e_flux[i_start];

    /* Interior cells */
    for (int i = i_start + 1; i < N - 1; i++) {
        rho_out[i] = rho_in[i] + dt / dx[i] * (rho_flux[i-1] - rho_flux[i]);
        double new_mom = mom[i] + dt / dx[i] * (mom_flux[i-1] - mom_flux[i]);
        u_out[i] = new_mom / rho_out[i];
        e_out[i] = e_in[i] + dt / dx[i] * (e_flux[i-1] - e_flux[i]);
    }

    /* Friction source term (implicit) */
    if (!g_no_source) {
        for (int i = i_start; i < N - 1; i++) {
            double denom_f = 1.0 + sqrt(1.0 + 2.0 * dt * fabs(u_out[i]) * lambda_fric / g_dia);
            u_out[i] = 2.0 * u_out[i] / denom_f;
        }
    }

    /* Heat transfer is now handled by update_wall_and_heat() in the main loop */

    /* Safety: cap velocity if internal energy went negative at outlet cell.
       Skip in tube mode (N-2 is just a regular cell near the right wall). */
    if (g_vacuum_ext <= 0.0 && u_out[N-2] > 0.0) {
        double eint_N2 = e_out[N-2] / rho_out[N-2] - 0.5 * u_out[N-2] * u_out[N-2];
        if (eint_N2 < 1000.0) {
            double e_total_sp = e_out[N-2] / rho_out[N-2];
            u_out[N-2] = sqrt(fmax(2.0 * (e_total_sp - 1000.0), 0.0));
            eint_N2 = e_total_sp - 0.5 * u_out[N-2] * u_out[N-2];
            if (eint_N2 < 1000.0) eint_N2 = 1000.0;
            e_out[N-2] = rho_out[N-2] * (eint_N2 + 0.5 * u_out[N-2] * u_out[N-2]);
        }
    }

    /* Copy ghost cell from input (main loop sets it before each HLLC call) */
    rho_out[N-1] = rho_in[N-1];
    u_out[N-1] = u_in[N-1];
    e_out[N-1] = e_in[N-1];
}

/* ========== Lax-Friedrich flux computation (forward Euler, no reconstruction) ========== */

static void compute_lf_step(
    const double *rho_in, const double *u_in, const double *e_in,
    const double *p_in, const double *a_in,
    const double *T_in, const double *Cv_in,
    int N, int i_start, double dt, double lambda_fric,
    /* outputs: */
    double *rho_out, double *u_out, double *e_out)
{
    double mom[MAX_CELLS];
    for (int i = i_start; i < N; i++) mom[i] = rho_in[i] * u_in[i];

    /* Lax-Friedrich fluxes at each face */
    int Nfaces = N - 1;

    /* Physical flux vectors at cell centres */
    /* F(U) = [rho*u, rho*u*u + p, u*(e + p)] */
    for (int f = i_start; f < Nfaces; f++) {
        int iL = f, iR = f + 1;

        double lambda_L = fabs(u_in[iL]) + a_in[iL];
        double lambda_R = fabs(u_in[iR]) + a_in[iR];
        double lam = fmax(lambda_L, lambda_R);

        /* F_rho */
        double FL_r = mom[iL];
        double FR_r = mom[iR];
        rho_flux[f] = 0.5 * (FL_r + FR_r) + 0.5 * lam * (rho_in[iL] - rho_in[iR]);

        /* F_mom */
        double FL_m = rho_in[iL] * u_in[iL] * u_in[iL] + p_in[iL];
        double FR_m = rho_in[iR] * u_in[iR] * u_in[iR] + p_in[iR];
        mom_flux[f] = 0.5 * (FL_m + FR_m) + 0.5 * lam * (mom[iL] - mom[iR]);

        /* F_energy */
        double FL_e = u_in[iL] * (e_in[iL] + p_in[iL]);
        double FR_e = u_in[iR] * (e_in[iR] + p_in[iR]);
        e_flux[f] = 0.5 * (FL_e + FR_e) + 0.5 * lam * (e_in[iL] - e_in[iR]);
    }

    /* Conservative update */
    /* Cell i_start (left wall BC) */
    rho_out[i_start] = rho_in[i_start] - dt / dx[i_start] * rho_flux[i_start];
    double mom1_val = rho_in[i_start] * u_in[i_start]
                    - dt / dx[i_start] * (mom_flux[i_start] - p_in[i_start]);
    u_out[i_start] = mom1_val / rho_out[i_start];
    e_out[i_start] = e_in[i_start] - dt / dx[i_start] * e_flux[i_start];

    /* Interior cells */
    for (int i = i_start + 1; i < N - 1; i++) {
        rho_out[i] = rho_in[i] + dt / dx[i] * (rho_flux[i-1] - rho_flux[i]);
        double new_mom = mom[i] + dt / dx[i] * (mom_flux[i-1] - mom_flux[i]);
        u_out[i] = new_mom / rho_out[i];
        e_out[i] = e_in[i] + dt / dx[i] * (e_flux[i-1] - e_flux[i]);
    }

    /* Friction source term (implicit) */
    if (!g_no_source) {
        for (int i = i_start; i < N - 1; i++) {
            double denom_f = 1.0 + sqrt(1.0 + 2.0 * dt * fabs(u_out[i]) * lambda_fric / g_dia);
            u_out[i] = 2.0 * u_out[i] / denom_f;
        }
    }

    /* Heat transfer is now handled by update_wall_and_heat() in the main loop */

    /* Safety: cap velocity if internal energy went negative at outlet cell.
       Skip in tube mode. */
    if (g_vacuum_ext <= 0.0 && u_out[N-2] > 0.0) {
        double eint_N2 = e_out[N-2] / rho_out[N-2] - 0.5 * u_out[N-2] * u_out[N-2];
        if (eint_N2 < 1000.0) {
            double e_total_sp = e_out[N-2] / rho_out[N-2];
            u_out[N-2] = sqrt(fmax(2.0 * (e_total_sp - 1000.0), 0.0));
            eint_N2 = e_total_sp - 0.5 * u_out[N-2] * u_out[N-2];
            if (eint_N2 < 1000.0) eint_N2 = 1000.0;
            e_out[N-2] = rho_out[N-2] * (eint_N2 + 0.5 * u_out[N-2] * u_out[N-2]);
        }
    }

    /* Copy ghost cell from input */
    rho_out[N-1] = rho_in[N-1];
    u_out[N-1] = u_in[N-1];
    e_out[N-1] = e_in[N-1];
}

/* ========== Recover thermodynamic state (Newton for T, then p, a) ========== */

static void recover_state(double *rho_in, double *u_in, double *e_in,
                          double *T_out, double *p_out, double *a_out,
                          double *Cv_out, int i_start, int N)
{
    for (int i = i_start; i < N; i++) {
        /* Guard: if conserved variables are NaN, reset to a safe state */
        if (rho_in[i] != rho_in[i] || e_in[i] != e_in[i] ||
            rho_in[i] <= 0.0) {
            rho_in[i] = g_pressure_R / (g_T0 * RGAS);
            u_in[i]   = 0.0;
            e_in[i]   = rho_in[i] * CV * g_T0;
            T_out[i]  = g_T0;
            Cv_out[i] = CV;
            p_out[i]  = g_pressure_R;
            a_out[i]  = 250.0;
            continue;
        }

        double e_int = (e_in[i] - 0.5 * rho_in[i] * u_in[i] * u_in[i]) / rho_in[i];

        if (g_ideal_gas) {
            /* Direct: no tables, no iteration */
            double Ti = e_int / CV;
            T_out[i] = Ti;
            Cv_out[i] = CV;
            p_out[i] = rho_in[i] * RGAS * Ti;
            a_out[i] = sqrt(GAMMA_I * p_out[i] / rho_in[i]);
        } else {
            double p_tab = p_out[i];
            if (p_tab < 1e5) p_tab = 1e5;

            double Ti = T_out[i];
            for (int iter = 0; iter < 3; iter++) {
                double cv = lookup_Cv(Ti, p_tab);
                Cv_out[i] = cv;
                Ti += (e_int - lookup_U(Ti, p_tab)) / cv;
            }

            if (Ti < 155.0 || Ti > 395.0 || Ti != Ti) {
                /* Outside table range or NaN — use ideal gas with CoolProp U offset */
                Ti = (e_int - g_U_offset) / CV;
                if (Ti < 50.0) Ti = 50.0;
                T_out[i] = Ti;
                Cv_out[i] = CV;
                p_out[i] = rho_in[i] * RGAS * Ti;
                a_out[i] = sqrt(GAMMA_I * p_out[i] / rho_in[i]);
            } else {
                T_out[i] = Ti;
                double cfac = compression_factor(Ti, p_out[i]);
                p_out[i] = Ti * rho_in[i] * RGAS * cfac;
                a_out[i] = thermo_soundspeed(Ti, p_out[i], rho_in[i], Cv_out[i]);
            }

            /* Final NaN guard on derived quantities */
            if (p_out[i] != p_out[i] || a_out[i] != a_out[i] || p_out[i] <= 0.0) {
                T_out[i] = fmax(e_int / CV, 50.0);
                Cv_out[i] = CV;
                p_out[i] = rho_in[i] * RGAS * T_out[i];
                a_out[i] = sqrt(GAMMA_I * p_out[i] / rho_in[i]);
            }
        }
    }
}

/* ========== Active front tracking ========== */

static void update_active_front(int N)
{
    /* Scan from left to find leftmost disturbed cell */
    int found = N - 1;
    for (int i = 0; i < N; i++) {
        if (fabs(p[i] - g_pressure_L) / g_pressure_L > 1e-6) {
            found = i;
            break;
        }
    }
    g_i_front = (found > 10) ? found - 10 : 0;
}

/* ========== Adaptive cell merging ========== */

static int merge_cells(int *pN, double merge_tol, double merge_dx_max)
{
    int N = *pN;
    int *merge_flag = (int *)calloc(N, sizeof(int));
    if (!merge_flag) { fprintf(stderr, "merge_cells: alloc failed\n"); return 0; }

    /* Mark pairs for merging (skip last 2 cells — outlet BC region) */
    for (int i = g_i_front; i <= N - 3; ) {
        double dp_val = fabs(p[i] - p[i+1]) / fmax(p[i], p[i+1]);
        double dT_val = fabs(T_arr[i] - T_arr[i+1]) / fmax(T_arr[i], T_arr[i+1]);
        double du_val = fabs(u[i] - u[i+1]) / fmax(fmax(fabs(u[i]), fabs(u[i+1])), 1.0);

        double max_diff = fmax(fmax(dp_val, dT_val), du_val);

        if (max_diff < merge_tol && dx[i] + dx[i+1] <= merge_dx_max) {
            merge_flag[i] = 1;
            i += 2;  /* skip the partner */
        } else {
            i += 1;
        }
    }

    /* In-place compaction pass */
    int j = 0;
    int merged = 0;
    for (int i = 0; i < N; ) {
        if (merge_flag[i] && i + 1 < N) {
            /* Conservative merge of cells i and i+1 */
            double m1 = rho[i] * dx[i], m2 = rho[i+1] * dx[i+1];
            double new_dx = dx[i] + dx[i+1];
            rho[j] = (m1 + m2) / new_dx;
            u[j]   = (m1 * u[i] + m2 * u[i+1]) / (m1 + m2);
            e[j]   = (e[i] * dx[i] + e[i+1] * dx[i+1]) / new_dx;
            dx[j]  = new_dx;
            T_arr[j] = 0.5 * (T_arr[i] + T_arr[i+1]);  /* initial guess for recovery */
            p[j]   = 0.5 * (p[i] + p[i+1]);              /* initial guess */
            /* Average wall temperatures of merged cells */
            for (int w = 0; w < N_WALL; w++)
                TW(j,w) = 0.5 * (TW(i,w) + TW(i+1,w));
            j++; i += 2; merged++;
        } else {
            if (j != i) {
                rho[j] = rho[i]; u[j] = u[i]; e[j] = e[i]; dx[j] = dx[i];
                T_arr[j] = T_arr[i]; p[j] = p[i]; a_arr[j] = a_arr[i]; Cv_real[j] = Cv_real[i];
                for (int w = 0; w < N_WALL; w++) TW(j,w) = TW(i,w);
            }
            j++; i += 1;
        }
    }
    *pN = j;
    N = j;

    /* Rebuild grid geometry */
    x_face[0] = -g_pipe_length;
    for (int i = 0; i < N; i++) {
        x_face[i+1] = x_face[i] + dx[i];
        xc[i] = 0.5 * (x_face[i] + x_face[i+1]);
    }
    x_face[N] = 0.0;  /* snap outlet */

    /* Full state recovery for merged cells */
    recover_state(rho, u, e, T_arr, p, a_arr, Cv_real, 0, N);

    /* Update g_i_front by rescanning */
    update_active_front(N);

    free(merge_flag);
    return merged;
}

/* ========== Outlet zone coarsening ========== */
/* One-time coarsening of fine outlet cells (dx_init) up to ~dx_max.
   Used after tube→outlet switch and after early output phase.
   Returns new cell count. */
static int coarsen_outlet_zone(int N)
{
    int N_old = N;

    /* Find the leftmost cell in the transition zone */
    int i_limit = 0;
    for (int i = N - 2; i >= 0; i--) {
        if (xc[i] < -g_dx_trans - 500.0) { i_limit = i + 1; break; }
    }

    /* Build coarsened blocks: work backward from outlet */
    int n_new = 0;
    int tmp_cap = N - i_limit + 10;
    double *r2  = (double*)malloc(tmp_cap * sizeof(double));
    double *u2  = (double*)malloc(tmp_cap * sizeof(double));
    double *e2  = (double*)malloc(tmp_cap * sizeof(double));
    double *d2  = (double*)malloc(tmp_cap * sizeof(double));
    double *T2  = (double*)malloc(tmp_cap * sizeof(double));
    double *p2  = (double*)malloc(tmp_cap * sizeof(double));
    double *tw2 = (double*)malloc(tmp_cap * N_WALL * sizeof(double));

    int i = N - 2;  /* last physical cell (ghost at N-1) */
    while (i >= i_limit) {
        double acc_m = 0, acc_mu = 0, acc_e = 0, acc_dx = 0;
        double T_lo = T_arr[i], T_hi = T_arr[i];
        double p_lo = p[i],     p_hi = p[i];
        double tw_sum[N_WALL];
        for (int w = 0; w < N_WALL; w++) tw_sum[w] = 0.0;

        while (i >= i_limit && acc_dx + dx[i] <= g_dx_max * 1.05) {
            double mk = rho[i] * dx[i];
            acc_m  += mk;
            acc_mu += mk * u[i];
            acc_e  += e[i] * dx[i];
            for (int w = 0; w < N_WALL; w++)
                tw_sum[w] += TW(i, w) * dx[i];
            acc_dx += dx[i];
            if (T_arr[i] < T_lo) T_lo = T_arr[i];
            if (T_arr[i] > T_hi) T_hi = T_arr[i];
            if (p[i] < p_lo) p_lo = p[i];
            if (p[i] > p_hi) p_hi = p[i];
            i--;
        }
        if (acc_dx < 1e-10) break;
        r2[n_new] = acc_m / acc_dx;
        u2[n_new] = acc_mu / acc_m;
        e2[n_new] = acc_e / acc_dx;
        d2[n_new] = acc_dx;
        T2[n_new] = 0.5 * (T_lo + T_hi);
        p2[n_new] = 0.5 * (p_lo + p_hi);
        for (int w = 0; w < N_WALL; w++)
            tw2[n_new * N_WALL + w] = tw_sum[w] / acc_dx;
        n_new++;
    }

    /* Write coarsened cells back: cells 0..i_limit-1 stay,
       then n_new coarsened cells (reversed), then ghost */
    int j = i_limit;
    for (int k = n_new - 1; k >= 0; k--) {
        rho[j] = r2[k]; u[j] = u2[k]; e[j] = e2[k]; dx[j] = d2[k];
        T_arr[j] = T2[k]; p[j] = p2[k];
        for (int w = 0; w < N_WALL; w++)
            TW(j, w) = tw2[k * N_WALL + w];
        j++;
    }
    N = j + 1;  /* +1 for ghost cell */
    free(r2); free(u2); free(e2); free(d2); free(T2); free(p2); free(tw2);

    /* Rebuild geometry */
    x_face[0] = -g_pipe_length;
    for (int ii = 0; ii < N; ii++) {
        x_face[ii+1] = x_face[ii] + dx[ii];
        xc[ii] = 0.5 * (x_face[ii] + x_face[ii+1]);
    }
    x_face[N] = 0.0;

    /* Recover state and set ghost */
    recover_state(rho, u, e, T_arr, p, a_arr, Cv_real, 0, N);
    set_outlet_ghost(rho, u, e, p, T_arr, a_arr, Cv_real, N);

    /* Update active front */
    update_active_front(N);

    /* Report */
    double dxmin = dx[0];
    for (int ii = 1; ii < N - 1; ii++)
        if (dx[ii] < dxmin) dxmin = dx[ii];
    printf("  [COARSEN] outlet zone regridded: %d -> %d cells, "
           "dx_min=%.1f m\n", N_old, N, dxmin);
    fflush(stdout);

    return N;
}

/* ========== Usage ========== */

static void print_usage(const char *prog)
{
    printf(
    "Usage: %s [OPTIONS]\n"
    "\n"
    "  1D compressible Euler solver for subsea natural gas pipeline rupture and\n"
    "  blowdown simulation with real-gas thermodynamics.\n"
    "\n"
    "Solver:\n"
    "  --lf                Use Lax-Friedrich / forward Euler (default: MUSCL-HLLC / SSP-RK2)\n"
    "  --ideal             Ideal gas mode (Z=1, bypass thermodynamic tables)\n"
    "\n"
    "Pipeline presets (set geometry, pressure, friction, grid, endtime):\n"
    "  --ns1ag             NS1 string A, German end  (Lubmin,    224 km, 164 bar)\n"
    "  --ns1ar             NS1 string A, Russian end (Portovaya, 999.8 km, 164 bar)\n"
    "  --ns1bg             NS1 string B, German end  (Lubmin,    217.7 km, 164 bar)\n"
    "  --ns1br             NS1 string B, Russian end (Portovaya, 1006.2 km, 164 bar)\n"
    "  --ns2ag             NS2 string A, German end  (Lubmin,    153.6 km, 103 bar)\n"
    "  --ns2ar             NS2 string A, Russian end (Portovaya, 1076.4 km, 103 bar)\n"
    "\n"
    "Physical parameters:\n"
    "  --pL <Pa>           Initial pipe pressure                [165e5]\n"
    "  --pR <Pa>           Ambient / outlet back-pressure       [7.1e5]\n"
    "  --T0 <K>            Initial gas & seawater temperature   [282]\n"
    "  --dia <m>           Pipe inner diameter                  [1.153]\n"
    "  --rough <m>         Wall roughness (Colebrook)           [1.7e-5]\n"
    "  --lambda <val>      Override Darcy friction factor        [computed]\n"
    "  --length <m>        Pipe length                          [4400]\n"
    "\n"
    "Source terms:\n"
    "  --nosrc             Disable friction and heat transfer\n"
    "  --noheat            Disable heat transfer only (keep friction)\n"
    "  --simple-wall       Simplified steady-state wall model\n"
    "  --nochoke           Disable choked outlet condition\n"
    "\n"
    "Grid:\n"
    "  --dx-init <m>       Cell size at outlet end              [1.0]\n"
    "  --dx-max <m>        Maximum cell size                    [0.5]\n"
    "  --dx-trans <m>      Transition length, fine to coarse    [1000]\n"
    "  --dx-left <m>       Fine cells at closed end (0=off)     [0]\n"
    "  --fac <val>         Geometric stretching factor          [1.0]\n"
    "  --merge             Enable adaptive cell merging\n"
    "  --no-merge          Disable adaptive cell merging\n"
    "  --merge-tol <f>     Merge tolerance (relative)           [0.01]\n"
    "  --merge-dx-max <m>  Max cell size after merging          [5000]\n"
    "\n"
    "Time and output:\n"
    "  --endtime <s>       Simulation end time                  [20]\n"
    "  --pstop <Pa>        Stop at this closed-end pressure     [disabled]\n"
    "  --nprofiles <n>     Number of profile snapshots          [100]\n"
    "  --prefix <str>      Output file prefix\n"
    "  --early-end <s>     Detailed early output phase end      [disabled]\n"
    "  --early-dt <s>      Early phase profile interval         [0.01]\n"
    "\n"
    "Tube mode (shock tube with vacuum extension):\n"
    "  --tube <m>          Vacuum extension length\n"
    "  --tube-start <s>    Use tube mode for this duration, then switch to outlet BC\n"
    "\n"
    "Output files:\n"
    "  timeseries_<prefix>[_lf].dat   Time history at outlet (13 columns)\n"
    "  profiles_<prefix>[_lf].dat     Spatial snapshots at regular intervals\n"
    "  early_<prefix>.dat             Detailed early profiles (if --early-end set)\n"
    "\n"
    "Examples:\n"
    "  %s --ns1ar                      Full NS1A Russian blowdown (HLLC)\n"
    "  %s --ns1ag --lf --pstop 8e5     NS1A German end, LF, stop at 8 bar\n"
    "  %s --length 1000 --pL 100e5     Custom 1 km pipe at 100 bar\n"
    "  %s --ns1ar --T0 279             Russian end, 279 K water temperature\n",
    prog, prog, prog, prog, prog);
}

/* ========== Main ========== */

int main(int argc, char **argv)
{
    int use_lf = 0;
    int no_source = 0;
    int ideal_gas = 0;
    const char *prefix = NULL;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
        else if (strcmp(argv[i], "--lf") == 0) use_lf = 1;
        else if (strcmp(argv[i], "--nosrc") == 0) no_source = 1;
        else if (strcmp(argv[i], "--noheat") == 0) g_no_heat = 1;
        else if (strcmp(argv[i], "--simple-wall") == 0) g_simple_wall = 1;
        else if (strcmp(argv[i], "--ideal") == 0) ideal_gas = 1;
        else if (strcmp(argv[i], "--nochoke") == 0) g_no_choke = 1;
        else if (strcmp(argv[i], "--merge") == 0) g_merge = 1;
        else if (strcmp(argv[i], "--ns1ag") == 0) {
            /* NS1 string A, German end (Lubmin): 224 km */
            g_pipe_length = 224000.0;  g_pressure_L = 164e5; g_pressure_R = 7.61e5;
            g_end_time = 50000.0; g_dx_init = 1.0; g_dx_max = 20.0;
            g_dx_trans = 1000.0; g_merge = 1; g_lambda_override = 0.007;
            if (!prefix) prefix = "ns1ag";
        }
        else if (strcmp(argv[i], "--ns1ar") == 0) {
            /* NS1 string A, Russian end (Portovaya): 999.8 km */
            g_pipe_length = 999800.0;  g_pressure_L = 164e5; g_pressure_R = 7.61e5;
            g_end_time = 450000.0; g_dx_init = 1.0; g_dx_max = 200.0;
            g_dx_trans = 1000.0; g_merge = 1; g_lambda_override = 0.007;
            if (!prefix) prefix = "ns1ar";
        }
        else if (strcmp(argv[i], "--ns1bg") == 0) {
            /* NS1 string B, German end (Lubmin): 217.7 km */
            g_pipe_length = 217700.0;  g_pressure_L = 164e5; g_pressure_R = 7.23e5;
            g_end_time = 50000.0; g_dx_init = 1.0; g_dx_max = 20.0;
            g_dx_trans = 1000.0; g_merge = 1; g_lambda_override = 0.007;
            if (!prefix) prefix = "ns1bg";
        }
        else if (strcmp(argv[i], "--ns1br") == 0) {
            /* NS1 string B, Russian end (Portovaya): 1006.2 km */
            g_pipe_length = 1006200.0;  g_pressure_L = 164e5; g_pressure_R = 7.23e5;
            g_end_time = 450000.0; g_dx_init = 1.0; g_dx_max = 200.0;
            g_dx_trans = 1000.0; g_merge = 1; g_lambda_override = 0.007;
            if (!prefix) prefix = "ns1br";
        }
        else if (strcmp(argv[i], "--ns2ag") == 0) {
            /* NS2 string A, German end (Lubmin): 153.6 km */
            g_pipe_length = 153600.0;  g_pressure_L = 103e5; g_pressure_R = 7.1e5;
            g_end_time = 50000.0; g_dx_init = 1.0; g_dx_max = 20.0;
            g_dx_trans = 1000.0; g_merge = 1;
            if (!prefix) prefix = "ns2ag";
        }
        else if (strcmp(argv[i], "--ns2ar") == 0) {
            /* NS2 string A, Russian end (Portovaya): 1076.4 km */
            g_pipe_length = 1076400.0;  g_pressure_L = 103e5; g_pressure_R = 7.1e5;
            g_end_time = 450000.0; g_dx_init = 1.0; g_dx_max = 200.0;
            g_dx_trans = 1000.0; g_merge = 1;
            if (!prefix) prefix = "ns2ar";
        }
        else if (strcmp(argv[i], "--length") == 0 && i+1 < argc)
            g_pipe_length = atof(argv[++i]);
        else if (strcmp(argv[i], "--endtime") == 0 && i+1 < argc)
            g_end_time = atof(argv[++i]);
        else if (strcmp(argv[i], "--fac") == 0 && i+1 < argc)
            g_fac = atof(argv[++i]);
        else if (strcmp(argv[i], "--dx-max") == 0 && i+1 < argc)
            g_dx_max = atof(argv[++i]);
        else if (strcmp(argv[i], "--dx-init") == 0 && i+1 < argc)
            g_dx_init = atof(argv[++i]);
        else if (strcmp(argv[i], "--dx-left") == 0 && i+1 < argc)
            g_dx_init_left = atof(argv[++i]);
        else if (strcmp(argv[i], "--prefix") == 0 && i+1 < argc)
            prefix = argv[++i];
        else if (strcmp(argv[i], "--nprofiles") == 0 && i+1 < argc)
            g_n_profiles = atoi(argv[++i]);
        else if (strcmp(argv[i], "--lambda") == 0 && i+1 < argc)
            g_lambda_override = atof(argv[++i]);
        else if (strcmp(argv[i], "--pL") == 0 && i+1 < argc)
            g_pressure_L = atof(argv[++i]);
        else if (strcmp(argv[i], "--pR") == 0 && i+1 < argc)
            g_pressure_R = atof(argv[++i]);
        else if (strcmp(argv[i], "--T0") == 0 && i+1 < argc)
            g_T0 = atof(argv[++i]);
        else if (strcmp(argv[i], "--dia") == 0 && i+1 < argc)
            g_dia = atof(argv[++i]);
        else if (strcmp(argv[i], "--rough") == 0 && i+1 < argc)
            g_rough = atof(argv[++i]);
        else if (strcmp(argv[i], "--no-merge") == 0)
            g_merge = 0;
        else if (strcmp(argv[i], "--merge-tol") == 0 && i+1 < argc)
            g_merge_tol = atof(argv[++i]);
        else if (strcmp(argv[i], "--merge-dx-max") == 0 && i+1 < argc)
            g_merge_dx_max = atof(argv[++i]);
        else if (strcmp(argv[i], "--dx-trans") == 0 && i+1 < argc)
            g_dx_trans = atof(argv[++i]);
        else if (strcmp(argv[i], "--early-end") == 0 && i+1 < argc)
            g_early_end = atof(argv[++i]);
        else if (strcmp(argv[i], "--early-dt") == 0 && i+1 < argc)
            g_early_dt = atof(argv[++i]);
        else if (strcmp(argv[i], "--pstop") == 0 && i+1 < argc)
            g_pstop = atof(argv[++i]);
        else if (strcmp(argv[i], "--tube") == 0 && i+1 < argc)
            g_vacuum_ext = atof(argv[++i]);
        else if (strcmp(argv[i], "--tube-start") == 0 && i+1 < argc) {
            g_tube_start = atof(argv[++i]);
            if (g_vacuum_ext <= 0.0) g_vacuum_ext = 2000.0; /* default 2km vacuum */
        }
    }
    g_no_source = no_source;
    g_ideal_gas = ideal_gas;

    /* Auto-enable early recording for tube-start runs to capture peak mass flow */
    if (g_tube_start > 0.0 && g_early_end <= 0.0) {
        g_early_end = g_tube_start + 3.0;  /* record 3s after tube→outlet switch */
        g_early_dt = 0.01;
    }

    /* Compute derived geometry from CLI parameters */
    g_area = g_dia * g_dia / 4.0 * PI_VAL;
    g_circum = g_dia * PI_VAL;
    /* Wall thermal model parameters printed below after grid setup */

    double lambda_fric = compute_lambda();
    if (g_lambda_override > 0.0) lambda_fric = g_lambda_override;
    printf("Friction: lambda = %.6f (Darcy)\n", lambda_fric);

    /* Compute CoolProp U offset for ideal-gas fallback.
       At low P, the gas is nearly ideal, so U_cp(T) ≈ CV*T + offset. */
    if (!g_ideal_gas) {
        g_U_offset = lookup_U(g_T0, 1e5) - CV * g_T0;
    }

    /* ===== Grid setup ===== */
    int n_faces;

    if (g_vacuum_ext > 0.0) {
        /* Tube mode: transition grid for pipe (-L to 0) + fine grid for vacuum (0 to +ext) */
        x_face[0] = -g_pipe_length;
        n_faces = 1;
        /* Pipe section: linear transition (same as non-tube mode) */
        while (x_face[n_faces - 1] < -1e-10) {
            double dist = -x_face[n_faces - 1];
            double dx_here;
            if (g_dx_trans > 0.0 && g_dx_max > g_dx_init) {
                if (dist > g_dx_trans) dx_here = g_dx_max;
                else dx_here = g_dx_init + (g_dx_max - g_dx_init) * dist / g_dx_trans;
                if (dx_here < g_dx_init) dx_here = g_dx_init;
            } else {
                dx_here = g_dx_init;
            }
            x_face[n_faces] = x_face[n_faces - 1] + dx_here;
            n_faces++;
            if (n_faces > MAX_CELLS - 5000) { fprintf(stderr, "Too many cells (pipe section)\n"); return 1; }
        }
        x_face[n_faces - 1] = 0.0;  /* snap outlet/diaphragm */
        /* If snapping created a tiny cell, merge it */
        if (n_faces > 2) {
            double last_dx = x_face[n_faces-1] - x_face[n_faces-2];
            if (last_dx < 0.5 * g_dx_init) {
                x_face[n_faces-2] = 0.0;
                n_faces--;
            }
        }
        /* Vacuum section: uniform at dx_init */
        int n_vac = (int)(g_vacuum_ext / g_dx_init + 0.5);
        for (int i = 1; i <= n_vac; i++) {
            x_face[n_faces] = i * g_dx_init;
            n_faces++;
        }
        x_face[n_faces - 1] = g_vacuum_ext;  /* snap end */
    } else if (g_dx_trans > 0.0 && !(g_fac > 1.0 && g_dx_init_left > 0.0)) {
        /* Linear transition grid: fine at outlet, linear ramp to dx_max */
        x_face[0] = -g_pipe_length;
        n_faces = 1;
        while (x_face[n_faces - 1] < -1e-10) {
            double dist = -x_face[n_faces - 1];  /* distance from outlet */
            double dx_here;
            if (dist > g_dx_trans)
                dx_here = g_dx_max;
            else
                dx_here = g_dx_init + (g_dx_max - g_dx_init) * dist / g_dx_trans;
            if (dx_here < g_dx_init) dx_here = g_dx_init;  /* floor */
            x_face[n_faces] = x_face[n_faces - 1] + dx_here;
            n_faces++;
            if (n_faces > MAX_CELLS) { fprintf(stderr, "Too many cells\n"); return 1; }
        }
        x_face[n_faces - 1] = 0.0;  /* snap outlet */
        /* If snapping created a tiny last cell, merge it with the previous one */
        if (n_faces > 2) {
            double last_dx = x_face[n_faces-1] - x_face[n_faces-2];
            if (last_dx < 0.5 * g_dx_init) {
                x_face[n_faces-2] = 0.0;
                n_faces--;
            }
        }
    } else if (g_dx_init_left > 0.0) {
        /* Two-sided stretching: fine at both ends, coarse in the middle (legacy). */
        n_faces = 1;
        x_face[0] = -g_pipe_length;
        double L = g_pipe_length;
        double slope = g_fac - 1.0;

        while (x_face[n_faces - 1] < -1e-10) {
            double x = x_face[n_faces - 1];
            double dist_outlet = -x;
            double dist_closed = x + L;

            double dx_right = g_dx_init + slope * dist_outlet;
            if (dx_right > g_dx_max) dx_right = g_dx_max;

            double dx_left = g_dx_init_left + slope * dist_closed;
            if (dx_left > g_dx_max) dx_left = g_dx_max;

            double dx_here = fmin(dx_right, dx_left);
            x_face[n_faces] = x + dx_here;
            n_faces++;
            if (n_faces > MAX_CELLS) { fprintf(stderr, "Too many cells\n"); return 1; }
        }
        x_face[n_faces - 1] = 0.0;  /* snap to exact outlet position */
    } else {
        /* One-sided stretching (original): fine at outlet only */
        n_faces = 1;
        x_face[0] = 0.0;
        double dx_c = -g_dx_init;

        while (x_face[n_faces - 1] > -g_pipe_length) {
            x_face[n_faces] = x_face[n_faces - 1] + dx_c;
            n_faces++;
            if (n_faces > MAX_CELLS) { fprintf(stderr, "Too many cells\n"); return 1; }
            if (-dx_c < g_dx_max) dx_c *= g_fac;
        }
        x_face[n_faces - 1] = -g_pipe_length;

        /* Reverse so x_face goes from -pipe_length to 0 */
        for (int i = 0; i < n_faces / 2; i++) {
            double tmp = x_face[i];
            x_face[i] = x_face[n_faces - 1 - i];
            x_face[n_faces - 1 - i] = tmp;
        }
    }

    int N = n_faces - 1;  /* number of cells */
    for (int i = 0; i < N; i++) {
        xc[i] = 0.5 * (x_face[i] + x_face[i+1]);
        dx[i] = x_face[i+1] - x_face[i];
    }

    printf("Mode: %s%s%s%s%s%s\n",
           use_lf ? "Lax-Friedrich (forward Euler)" : "MUSCL-HLLC (SSP-RK2)",
           ideal_gas ? " [IDEAL GAS]" : "",
           no_source ? " [NO SOURCE]" : "",
           g_no_choke ? " [NO CHOKE]" : "",
           g_merge ? " [MERGE]" : "",
           g_vacuum_ext > 0.0 ? " [TUBE]" : "");
    printf("Pipe: %.1f km, endTime=%.0f s, dx_init=%.1f m, dx_max=%.1f m",
           g_pipe_length / 1000.0, g_end_time, g_dx_init, g_dx_max);
    if (g_dx_trans > 0.0 && !(g_fac > 1.0 && g_dx_init_left > 0.0))
        printf(", dx_trans=%.0f m", g_dx_trans);
    else
        printf(", fac=%.4f", g_fac);
    printf("\n");
    printf("Conditions: pL=%.1f bar, pR=%.2f bar, T0=%.1f K, dia=%.3f m\n",
           g_pressure_L / 1e5, g_pressure_R / 1e5, g_T0, g_dia);
    double dx_min_actual = dx[0], dx_max_actual = dx[0];
    for (int i = 1; i < N; i++) {
        if (dx[i] < dx_min_actual) dx_min_actual = dx[i];
        if (dx[i] > dx_max_actual) dx_max_actual = dx[i];
    }
    printf("Grid: %d cells, %d faces, dx_min=%.2f m, dx_max=%.2f m\n",
           N, n_faces, dx_min_actual, dx_max_actual);
    if (g_merge)
        printf("Merge: tol=%.4f, dx_max=%.0f m\n", g_merge_tol, g_merge_dx_max);

    /* Initialize active front for large grids (skip in tube mode) */
    if (g_vacuum_ext <= 0.0 && N > 200) {
        for (int i = N - 1; i >= 0; i--) {
            if (-xc[i] > g_dx_trans + 500.0) { g_i_front = i; break; }
        }
        if (g_i_front > 20) g_i_front -= 20;  /* safety margin */
        else g_i_front = 0;
        printf("Active front: i_front=%d (x=%.0f m)\n", g_i_front, xc[g_i_front]);
    }

    /* Find diaphragm cell in tube mode */
    if (g_vacuum_ext > 0.0) {
        for (int i = 0; i < N; i++) {
            if (xc[i] >= 0.0) { g_i_diaphragm = i; break; }
        }
        printf("Tube mode: %.0f m pressurized + %.0f m vacuum, diaphragm at cell %d (x=%.1f m)\n",
               g_pipe_length, g_vacuum_ext, g_i_diaphragm, xc[g_i_diaphragm]);
    }
    fflush(stdout);

    /* ===== Initial conditions ===== */
    double cfacR = compression_factor(g_T0, g_pressure_R);
    double cfacL = compression_factor(g_T0, g_pressure_L);
    double densityR = g_pressure_R / (g_T0 * RGAS * cfacR);
    double densityL = g_pressure_L / (g_T0 * RGAS * cfacL);

    if (g_vacuum_ext > 0.0) {
        /* Tube mode: left side at pL, right side at pR */
        for (int i = 0; i < N; i++) {
            if (xc[i] < 0.0) {
                rho[i] = densityL;
                p[i]   = g_pressure_L;
            } else {
                rho[i] = densityR;
                p[i]   = g_pressure_R;
            }
            u[i] = 0.0;
            T_arr[i] = g_T0;
            e[i] = rho[i] * lookup_U(g_T0, p[i]);
        }
    } else {
        /* Normal mode: all interior cells at pL, ghost cell at pR */
        for (int i = 0; i < N; i++) {
            if (i < N - 1) {
                rho[i] = densityL;
                p[i]   = g_pressure_L;
            } else {
                rho[i] = densityR;
                p[i]   = g_pressure_R;
            }
            u[i] = 0.0;
            T_arr[i] = g_T0;
            e[i] = rho[i] * lookup_U(g_T0, p[i]);
        }
    }

    /* Initialize wall temperatures to T0 */
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N_WALL; j++)
            TW(i,j) = g_T0;

    if (g_simple_wall)
        printf("Wall model: simplified (steady-state, T_wall=%.0f K, R=%.4f m2K/W)\n",
               g_T0, T_STEEL_WALL / K_STEEL + T_CONC_WALL / K_CONC);
    else
        printf("Wall model: %d steel nodes (%.1f mm each) + %d concrete nodes (%.1f mm each)\n",
               N_WALL_STEEL, T_STEEL_WALL / N_WALL_STEEL * 1000,
               N_WALL_CONC, T_CONC_WALL / N_WALL_CONC * 1000);

    /* ===== Output setup ===== */
    double dt_tseries = g_end_time / N_TSERIES;
    double dt_profiles = g_end_time / g_n_profiles;
    double next_tseries = 0.0;
    double next_profiles = 0.0;

    char ts_file[256], pr_file[256], early_file[256];
    if (prefix) {
        snprintf(ts_file, sizeof(ts_file), "timeseries_%s%s.dat",
                 prefix, use_lf ? "_lf" : "");
        snprintf(pr_file, sizeof(pr_file), "profiles_%s%s.dat",
                 prefix, use_lf ? "_lf" : "");
        snprintf(early_file, sizeof(early_file), "early_%s%s.dat",
                 prefix, use_lf ? "_lf" : "");
    } else {
        snprintf(ts_file, sizeof(ts_file), "%s",
                 use_lf ? "timeseries_lf.dat" : "timeseries.dat");
        snprintf(pr_file, sizeof(pr_file), "%s",
                 use_lf ? "profiles_lf.dat" : "profiles.dat");
        snprintf(early_file, sizeof(early_file), "%s",
                 use_lf ? "early_lf.dat" : "early.dat");
    }
    FILE *fp_ts = fopen(ts_file, "w");
    FILE *fp_pr = fopen(pr_file, "w");
    if (!fp_ts || !fp_pr) { fprintf(stderr, "Cannot open output files\n"); return 1; }
    printf("Output: %s, %s\n", ts_file, pr_file);

    FILE *fp_early = NULL;
    double next_early = 0.0;
    if (g_early_end > 0.0) {
        fp_early = fopen(early_file, "w");
        if (!fp_early) { fprintf(stderr, "Cannot open early output file\n"); return 1; }
        printf("Early output: %s (until t=%.2f s, dt=%.4f s)\n",
               early_file, g_early_end, g_early_dt);
        fprintf(fp_early, "# Early detailed profiles (x > -1000 m)\n");
    }

    fflush(stdout);
    if (g_vacuum_ext > 0.0)
        fprintf(fp_ts, "# t massflow[kg/s] u_diaph rho_diaph p_diaph T_diaph totmass p_base p_mid a_diaph Mach_diaph poutlet_choked choked_active\n");
    else
        fprintf(fp_ts, "# t massflow[kg/s] u_out rho_out p_out T_out totmass p_base p_mid a_out Mach_out poutlet_choked choked_active\n");

    /* ===== Time integration ===== */
    double sim_time = 0.0;
    int nstep = 0;
    struct timespec ts_start, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    /* Initialise a_arr for first CFL */
    for (int i = 0; i < N; i++) {
        a_arr[i] = thermo_soundspeed(T_arr[i], p[i], rho[i], CV);
    }

    /* Initial state recovery (before first CFL and output) */
    recover_state(rho, u, e, T_arr, p, a_arr, Cv_real, 0, N);

    /* Set initial outlet ghost cell (choked or subsonic) */
    set_outlet_ghost(rho, u, e, p, T_arr, a_arr, Cv_real, N);

    /* Merge schedule */
    double next_merge = 0.1;  /* first merge at t=0.1s */

    while (sim_time <= g_end_time) {
        nstep++;

        /* Advance active front to stay ahead of wave propagation.
           The wave moves at most CFL_NUM cells per step (~0.6), so
           retreating by 2 cells per step is always safe. */
        if (g_i_front > 0) {
            g_i_front = (g_i_front > 2) ? g_i_front - 2 : 0;
        }

        /* CFL timestep (active cells only) */
        double dt = 1e30;
        for (int i = g_i_front; i < N; i++) {
            double local_dt = CFL_NUM / (a_arr[i] + fabs(u[i])) * dx[i];
            if (local_dt < dt) dt = local_dt;
        }
        sim_time += dt;
        g_sim_time = sim_time;

        /* --- Tube-to-outlet BC switch --- */
        if (g_tube_start > 0.0 && g_vacuum_ext > 0.0 && sim_time >= g_tube_start) {
            /* Trim vacuum cells: keep only cells with x_face < ~0 */
            int N_new = g_i_diaphragm;  /* last pipe cell is i_diaphragm - 1 */
            if (N_new < 10) N_new = 10;
            /* The new ghost cell (N_new - 1) is at x ≈ 0 (outlet) */
            /* Copy the physical state of the old diaphragm cell into the ghost */
            N_new++;  /* +1 for ghost cell */
            N = N_new;
            /* Recompute dx for the trimmed grid */
            for (int i = 0; i < N; i++)
                dx[i] = x_face[i+1] - x_face[i];
            /* Switch off tube mode */
            g_vacuum_ext = 0.0;
            g_i_diaphragm = -1;
            g_tube_start = 0.0;
            /* Recover state and set outlet ghost */
            recover_state(rho, u, e, T_arr, p, a_arr, Cv_real, 0, N);
            set_outlet_ghost(rho, u, e, p, T_arr, a_arr, Cv_real, N);
            /* Initialize active front for the (now large) grid */
            if (N > 200) {
                for (int i = N - 1; i >= 0; i--) {
                    if (-xc[i] > g_dx_trans + 500.0) { g_i_front = i; break; }
                }
                if (g_i_front > 20) g_i_front -= 20;
                else g_i_front = 0;
            }
            /* Enable merge for the long run */
            if (!g_merge) g_merge = 1;
            printf("  [TUBE->OUTLET t=%.2f] Switched to outlet BC: N=%d, i_front=%d\n",
                   sim_time, N, g_i_front);
            fflush(stdout);

            /* Immediately coarsen fine outlet cells to speed up blowdown */
            N = coarsen_outlet_zone(N);
        }

        if (use_lf) {
            /* ===== Lax-Friedrich: single forward Euler step ===== */
            compute_lf_step(rho, u, e, p, a_arr, T_arr, Cv_real,
                            N, g_i_front, dt, lambda_fric,
                            rho1, u1, e1);
            /* Copy interior cells only (ghost is set after recover_state) */
            for (int i = g_i_front; i < N - 1; i++) {
                rho[i] = rho1[i];
                u[i]   = u1[i];
                e[i]   = e1[i];
            }
        } else {
            /* ===== HLLC with SSP-RK2 ===== */

            /* Save stage-0 state */
            memcpy(rho_0 + g_i_front, rho + g_i_front, (N - g_i_front) * sizeof(double));
            for (int i = g_i_front; i < N; i++) mom_0[i] = rho[i] * u[i];
            memcpy(e_0 + g_i_front, e + g_i_front, (N - g_i_front) * sizeof(double));

            /* Stage 1 */
            compute_hllc_flux(rho, u, e, p, T_arr, Cv_real,
                              N, g_i_front, dt, lambda_fric, 0,
                              rho1, u1, e1);

            /* Recover state after stage 1 */
            memcpy(T1 + g_i_front, T_arr + g_i_front, (N - g_i_front) * sizeof(double));
            memcpy(p1_arr + g_i_front, p + g_i_front, (N - g_i_front) * sizeof(double));
            recover_state(rho1, u1, e1, T1, p1_arr, a1_arr, Cv_real1, g_i_front, N);

            /* Set ghost cell for stage 2 using recovered stage-1 state */
            set_outlet_ghost(rho1, u1, e1, p1_arr, T1, a1_arr, Cv_real1, N);

            /* Stage 2 */
            compute_hllc_flux(rho1, u1, e1, p1_arr, T1, Cv_real1,
                              N, g_i_front, dt, lambda_fric, 1,
                              rho2, u2_arr, e2);

            /* SSP-RK2 combination: 0.5*stage0 + 0.5*stage2 */
            for (int i = g_i_front; i < N - 1; i++) {
                rho[i] = 0.5 * rho_0[i] + 0.5 * rho2[i];
                e[i]   = 0.5 * e_0[i]   + 0.5 * e2[i];
                double new_mom = 0.5 * mom_0[i] + 0.5 * (rho2[i] * u2_arr[i]);
                u[i] = new_mom / rho[i];
            }

        }

        /* Wall conduction + heat exchange (uses T_arr from previous step as BC) */
        update_wall_and_heat(N, g_i_front, e, T_arr, Cv_real, rho, u, dt);

        /* Recover thermodynamic state after update (for CFL, output, next step) */
        recover_state(rho, u, e, T_arr, p, a_arr, Cv_real, g_i_front, N);

        /* Clamp near-outlet cells to subsonic using exact recovered sound speed.
           For ideal gas this never fires; for real gas the approximate
           isentropic ghost relations can push M slightly > 1 in a band
           of cells near the outlet.  Apply to last 20 physical cells.
           Skip in tube mode (no outlet BC artifact). */
        if (!g_ideal_gas && g_vacuum_ext <= 0.0) {
            int i_clamp_start = (N - 2 > 20) ? N - 22 : 0;
            int clamped = 0;
            for (int i = i_clamp_start; i < N - 1; i++) {
                if (u[i] > a_arr[i]) {
                    u[i] = 0.999 * a_arr[i];
                    double eint = e[i] / rho[i] - 0.5 * u[i] * u[i];
                    if (eint < 1000.0) eint = 1000.0;
                    e[i] = rho[i] * (eint + 0.5 * u[i] * u[i]);
                    clamped = 1;
                }
            }
            if (clamped)
                recover_state(rho, u, e, T_arr, p, a_arr, Cv_real,
                              i_clamp_start, N - 1);
        }

        /* Set ghost cell for next timestep */
        set_outlet_ghost(rho, u, e, p, T_arr, a_arr, Cv_real, N);

        /* ===== Adaptive merge check ===== */
        if (g_merge && g_vacuum_ext <= 0.0 && sim_time >= next_merge) {
            update_active_front(N);

            double mass_before = 0.0;
            for (int i = 0; i < N; i++) mass_before += rho[i] * dx[i];

            int merged = merge_cells(&N, g_merge_tol, g_merge_dx_max);

            if (merged > 0) {
                double mass_after = 0.0;
                for (int i = 0; i < N; i++) mass_after += rho[i] * dx[i];

                double dx_min_m = dx[0], dx_max_m = dx[0];
                for (int i = 1; i < N; i++) {
                    if (dx[i] < dx_min_m) dx_min_m = dx[i];
                    if (dx[i] > dx_max_m) dx_max_m = dx[i];
                }
                printf("  [MERGE t=%.2f] %d pairs merged -> %d cells, "
                       "dx_min=%.1f, dx_max=%.1f, i_front=%d, mass_err=%.2e\n",
                       sim_time, merged, N, dx_min_m, dx_max_m, g_i_front,
                       (mass_after - mass_before) / mass_before);
                fflush(stdout);
            }

            next_merge = sim_time * 2.0;
            if (next_merge < sim_time + 1.0) next_merge = sim_time + 1.0;
        }

        /* Periodic active front rescan (independent of merge) */
        if (g_i_front > 0 && nstep % 10000 == 0) {
            update_active_front(N);
        }

        /* ===== Output ===== */
        if (sim_time >= next_tseries) {
            /* During early phase, write timeseries at early_dt for peak resolution */
            if (fp_early && g_early_end > 0.0 && sim_time < g_early_end)
                next_tseries += g_early_dt;
            else
                next_tseries += dt_tseries;
            double totmass = 0.0;
            for (int i = 0; i < N; i++) totmass += rho[i] * dx[i] * g_area;
            /* In tube mode, report values at the diaphragm (x=0).
               In normal mode, report outlet ghost cell (choked/subsonic BC). */
            int ir = (g_vacuum_ext > 0.0 && g_i_diaphragm >= 0) ? g_i_diaphragm : N-1;
            double mf = u[ir] * rho[ir] * g_area;
            int i_mid = N / 2;
            double a_out = a_arr[ir];
            double mach_out = (a_out > 0.0) ? u[ir] / a_out : 0.0;
            fprintf(fp_ts, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %d\n",
                    sim_time, mf, u[ir], rho[ir], p[ir], T_arr[ir], totmass, p[0], p[i_mid], a_out, mach_out,
                    g_poutlet_choked, g_choked_active);
        }

        if (sim_time >= next_profiles) {
            next_profiles += dt_profiles;
            fprintf(fp_pr, "# t = %.8e\n", sim_time);
            for (int i = 0; i < N; i++) {
                double e_int = (e[i] - 0.5 * rho[i] * u[i] * u[i]) / rho[i];
                fprintf(fp_pr, "%.6e %.8e %.8e %.8e %.4f %.8e\n",
                        xc[i], rho[i], u[i], p[i], T_arr[i], e_int);
            }
            fprintf(fp_pr, "\n");

            clock_gettime(CLOCK_MONOTONIC, &ts_now);
            double elapsed = (ts_now.tv_sec - ts_start.tv_sec)
                           + (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
            printf("t = %.4f / %.1f  steps = %d  cells = %d  steps/s = %.1f\n",
                   sim_time, g_end_time, nstep, N, nstep / elapsed);
            fflush(stdout);
        }

        /* Early detailed output (near outlet only) */
        if (fp_early && sim_time < g_early_end && sim_time >= next_early) {
            next_early += g_early_dt;
            fprintf(fp_early, "# t = %.8e\n", sim_time);
            for (int i = 0; i < N; i++) {
                if (xc[i] > -1000.0) {
                    double e_int = (e[i] - 0.5 * rho[i] * u[i] * u[i]) / rho[i];
                    fprintf(fp_early, "%.6e %.8e %.8e %.8e %.4f %.8e\n",
                            xc[i], rho[i], u[i], p[i], T_arr[i], e_int);
                }
            }
            fprintf(fp_early, "\n");
        }

        /* Close early output in tube mode (no coarsening needed) */
        if (fp_early && sim_time >= g_early_end && g_vacuum_ext > 0.0) {
            fclose(fp_early);
            fp_early = NULL;
        }

        /* One-time outlet coarsening after early phase completes.
           Skip in tube mode (uniform grid, no transition zone). */
        if (fp_early && sim_time >= g_early_end && g_vacuum_ext <= 0.0) {
            fclose(fp_early);
            fp_early = NULL;  /* prevents re-entry */
            N = coarsen_outlet_zone(N);
        }

        /* Stop when closed-end pressure drops below threshold */
        if (g_pstop > 0.0 && p[0] < g_pstop) {
            printf("Stopping: p_base = %.2f bar < pstop = %.2f bar\n",
                   p[0] / 1e5, g_pstop / 1e5);
            break;
        }
    }

    fclose(fp_ts);
    fclose(fp_pr);
    if (fp_early) fclose(fp_early);

    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    double elapsed = (ts_now.tv_sec - ts_start.tv_sec)
                   + (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
    printf("\nDone: %d steps in %.2f seconds (%.0f steps/s)\n", nstep, elapsed, nstep / elapsed);

    /* Mass conservation check */
    double totmass = 0.0;
    for (int i = 0; i < N; i++) totmass += rho[i] * dx[i] * g_area;
    printf("Final total mass: %.6f kg\n", totmass);
    printf("Final cells: %d\n", N);

    return 0;
}
