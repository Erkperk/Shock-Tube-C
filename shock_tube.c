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
#define CV    1622.0
#define CP    2125.0
#define RGAS  (8314.46 / 16.71)
#define GAMMA (CP / CV)

/* Pipe geometry */
#define DIA       1.153
#define ROUGH     (0.017e-3)
#define PI_VAL    3.14159265358979323846
#define AREA      (DIA * DIA / 4.0 * PI_VAL)
#define CIRCUM    (DIA * PI_VAL)
#define K_CONCRETE 4.0
#define THICK     (4.3 * 2.54 * 0.01)
#define WATT_PER_K_PER_M  (K_CONCRETE / THICK * CIRCUM)
#define WATT_PER_K_PER_M3 (WATT_PER_K_PER_M / AREA)

/* Numerics */
#define CFL_NUM    0.6
#define END_TIME   20.0
#define DX_INIT    1.0
#define DX_MAX     0.5
#define FAC        1.0

/* Pipeline */
#define PIPE_LENGTH  4400.0
#define PRESSURE_R   (7.1e5)
#define PRESSURE_L   (165e5)
#define T0           282.0

/* Output */
#define N_TSERIES  5000
#define N_PROFILES 100

/* Maximum number of cells */
#define MAX_CELLS 100000

/* Global flag to disable source terms for diagnostics */
static int g_no_source = 0;
/* Global flag for ideal gas mode (bypass tables, Z=1) */
static int g_ideal_gas = 0;

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
    return table_lookup(Cv_table, T, P);
}
static inline double lookup_a(double T, double P)   {
    if (g_ideal_gas) { (void)P; return sqrt(GAMMA * RGAS * T); }
    return table_lookup(a_table, T, P);
}

/* ========== Compression factor (analytical, 19 coefficients) ========== */

static inline double compression_factor(double T, double P)
{
    if (g_ideal_gas) { (void)T; (void)P; return 1.0; }
    if (T < 150.0) T = 150.0;

    const double gammag = 0.577;
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
    if (g_ideal_gas) return sqrt(GAMMA * P / rho_val);

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

/* ========== MC limiter ========== */

static inline double mc_limiter(double a, double b)
{
    if (a * b <= 0.0) return 0.0;
    double sa = (a > 0.0) ? 1.0 : -1.0;
    double aa = fabs(a), ab = fabs(b);
    double m = 2.0 * aa;
    if (2.0 * ab < m) m = 2.0 * ab;
    double c = 0.5 * fabs(a + b);
    if (c < m) m = c;
    return sa * m;
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
    double x = -2.0 * log10(ROUGH / (3.7 * DIA));
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

/* ========== HLLC flux computation ========== */

static void compute_hllc_flux(
    const double *rho_in, const double *u_in, const double *e_in,
    const double *p_in, const double *a_in,
    const double *T_in, const double *Cv_in,
    int N, double dt, double lambda_fric,
    int is_stage2 __attribute__((unused)),
    /* outputs: */
    double *rho_out, double *u_out, double *e_out)
{
    double mom[MAX_CELLS];
    for (int i = 0; i < N; i++) mom[i] = rho_in[i] * u_in[i];

    /* MUSCL reconstruction with van Leer limiter on CONSERVED variables */
    s_rho[0] = 0.0; s_mom[0] = 0.0; s_e[0] = 0.0;
    s_rho[N-1] = 0.0; s_mom[N-1] = 0.0; s_e[N-1] = 0.0;
    for (int i = 1; i < N - 1; i++) {
        double dxL = xc[i] - xc[i-1];
        double dxR = xc[i+1] - xc[i];
        s_rho[i] = vanleer_limiter((rho_in[i] - rho_in[i-1]) / dxL, (rho_in[i+1] - rho_in[i]) / dxR);
        s_mom[i] = vanleer_limiter((mom[i] - mom[i-1]) / dxL, (mom[i+1] - mom[i]) / dxR);
        s_e[i]   = vanleer_limiter((e_in[i] - e_in[i-1]) / dxL, (e_in[i+1] - e_in[i]) / dxR);
    }

    /* Compute HLLC fluxes at each face */
    int Nfaces = N - 1;
    for (int f = 0; f < Nfaces; f++) {
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
            /* Direct: p = (gamma-1) * rho * e_int, a = sqrt(gamma*p/rho) */
            pL = (GAMMA - 1.0) * rhoL * e_int_L;
            pR = (GAMMA - 1.0) * rhoR * e_int_R;
            aL = sqrt(GAMMA * pL / rhoL);
            aR = sqrt(GAMMA * pR / rhoR);
        } else {
            /* Iterated face recovery: Newton for T + update P each iteration */
            double TfL = T_in[iL];
            double PfL = fmax(p_in[iL], 1e5);
            for (int it = 0; it < 3; it++) {
                TfL += (e_int_L - lookup_U(TfL, PfL)) / lookup_Cv(TfL, PfL);
                if (TfL < 150.0) TfL = 150.0;
                if (TfL > 400.0) TfL = 400.0;
                PfL = rhoL * compression_factor(TfL, PfL) * RGAS * TfL;
            }
            pL = PfL;
            aL = thermo_soundspeed(TfL, PfL, rhoL, lookup_Cv(TfL, PfL));

            double TfR = T_in[iR];
            double PfR = fmax(p_in[iR], 1e5);
            for (int it = 0; it < 3; it++) {
                TfR += (e_int_R - lookup_U(TfR, PfR)) / lookup_Cv(TfR, PfR);
                if (TfR < 150.0) TfR = 150.0;
                if (TfR > 400.0) TfR = 400.0;
                PfR = rhoR * compression_factor(TfR, PfR) * RGAS * TfR;
            }
            pR = PfR;
            aR = thermo_soundspeed(TfR, PfR, rhoR, lookup_Cv(TfR, PfR));
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

    /* Choked outlet BC */
    {
        double a_ideal = sqrt(GAMMA * p_in[N-2] / rho_in[N-2]);

        if (u_in[N-2] < a_ideal && u_in[N-2] > 0.0) {
            double vc = u_in[N-2] / a_ideal;
            double xd = vc * (1.0 + GAMMA) / (2.0 + vc * (GAMMA - 1.0));
            double uoutlet = u_in[N-2] / xd;
            double c0 = 0.5 * (1.0 + GAMMA) * uoutlet;
            double p0 = p_in[N-2] * pow(c0 / a_ideal, 2.0 * GAMMA / (GAMMA - 1.0));
            double poutlet = p0 * pow(2.0 / (1.0 + GAMMA), 2.0 * GAMMA / (GAMMA - 1.0));
            double rhooutlet = GAMMA * poutlet / (uoutlet * uoutlet);
            double Toutlet = poutlet / (rhooutlet * RGAS);
            double T_clamped = fmax(Toutlet, 150.0);
            double P_clamped = fmax(poutlet, 1e5);
            double eoutlet = 0.5 * rhooutlet * uoutlet * uoutlet
                           + rhooutlet * lookup_U(T_clamped, P_clamped);

            if (poutlet > PRESSURE_R) {
                rho_flux[Nfaces-1] = rhooutlet * uoutlet;
                e_flux[Nfaces-1]   = uoutlet * (eoutlet + poutlet);
                mom_flux[Nfaces-1] = rhooutlet * uoutlet * uoutlet + poutlet;
            }
        } else if (u_in[N-2] == 0.0) {
            double uoutlet = 2.0 / (1.0 + GAMMA) * a_ideal;
            double poutlet = p_in[N-2] * pow(2.0 / (1.0 + GAMMA), 2.0 * GAMMA / (GAMMA - 1.0));
            double rhooutlet = GAMMA * poutlet / (uoutlet * uoutlet);
            double Toutlet = poutlet / (rhooutlet * RGAS);
            double T_clamped = fmax(Toutlet, 150.0);
            double P_clamped = fmax(poutlet, 1e5);
            double eoutlet = 0.5 * rhooutlet * uoutlet * uoutlet
                           + rhooutlet * lookup_U(T_clamped, P_clamped);

            rho_flux[Nfaces-1] = rhooutlet * uoutlet;
            e_flux[Nfaces-1]   = uoutlet * (eoutlet + poutlet);
            mom_flux[Nfaces-1] = rhooutlet * uoutlet * uoutlet + poutlet;
        }
    }

    /* Conservative update */
    /* Cell 0 (left boundary: wall/reflective — pressure contributes) */
    rho_out[0] = rho_in[0] - dt / dx[0] * rho_flux[0];
    double mom1_val = rho_in[0] * u_in[0] - dt / dx[0] * (mom_flux[0] - p_in[0]);
    u_out[0] = mom1_val / rho_out[0];
    e_out[0] = e_in[0] - dt / dx[0] * e_flux[0];

    /* Interior cells */
    for (int i = 1; i < N - 1; i++) {
        rho_out[i] = rho_in[i] + dt / dx[i] * (rho_flux[i-1] - rho_flux[i]);
        double new_mom = mom[i] + dt / dx[i] * (mom_flux[i-1] - mom_flux[i]);
        u_out[i] = new_mom / rho_out[i];
        e_out[i] = e_in[i] + dt / dx[i] * (e_flux[i-1] - e_flux[i]);
    }

    /* Friction source term (implicit) */
    if (!g_no_source) {
        for (int i = 0; i < N; i++) {
            double denom_f = 1.0 + sqrt(1.0 + 2.0 * dt * fabs(u_out[i]) * lambda_fric / DIA);
            u_out[i] = 2.0 * u_out[i] / denom_f;
        }
    }

    /* Heat transfer source term */
    if (!g_no_source) {
        for (int i = 0; i < N; i++) {
            double tdiff = T0 - T_in[i];
            double deltaT = tdiff * (1.0 - exp(-dt * WATT_PER_K_PER_M3 / (rho_in[i] * Cv_in[i])));
            double deltaE = deltaT * rho_in[i] * Cv_in[i];
            e_out[i] += deltaE;
        }
    }

    /* Right boundary */
    double cfac_bc = compression_factor(T0, PRESSURE_R);
    u_out[N-1] = u_out[N-2];
    rho_out[N-1] = PRESSURE_R / (cfac_bc * RGAS * T0);
    e_out[N-1] = rho_out[N-1] * lookup_U(T0, PRESSURE_R)
               + 0.5 * rho_out[N-1] * u_out[N-1] * u_out[N-1];
}

/* ========== Lax-Friedrich flux computation (forward Euler, no reconstruction) ========== */

static void compute_lf_step(
    const double *rho_in, const double *u_in, const double *e_in,
    const double *p_in, const double *a_in,
    const double *T_in, const double *Cv_in,
    int N, double dt, double lambda_fric,
    /* outputs: */
    double *rho_out, double *u_out, double *e_out)
{
    double mom[MAX_CELLS];
    for (int i = 0; i < N; i++) mom[i] = rho_in[i] * u_in[i];

    /* Lax-Friedrich fluxes at each face */
    int Nfaces = N - 1;

    /* Physical flux vectors at cell centres */
    /* F(U) = [rho*u, rho*u*u + p, u*(e + p)] */
    for (int f = 0; f < Nfaces; f++) {
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

    /* Choked outlet BC (same as HLLC version) */
    {
        double a_ideal = sqrt(GAMMA * p_in[N-2] / rho_in[N-2]);

        if (u_in[N-2] < a_ideal && u_in[N-2] > 0.0) {
            double vc = u_in[N-2] / a_ideal;
            double xd = vc * (1.0 + GAMMA) / (2.0 + vc * (GAMMA - 1.0));
            double uoutlet = u_in[N-2] / xd;
            double c0 = 0.5 * (1.0 + GAMMA) * uoutlet;
            double p0 = p_in[N-2] * pow(c0 / a_ideal, 2.0 * GAMMA / (GAMMA - 1.0));
            double poutlet = p0 * pow(2.0 / (1.0 + GAMMA), 2.0 * GAMMA / (GAMMA - 1.0));
            double rhooutlet = GAMMA * poutlet / (uoutlet * uoutlet);
            double Toutlet = poutlet / (rhooutlet * RGAS);
            double T_clamped = fmax(Toutlet, 150.0);
            double P_clamped = fmax(poutlet, 1e5);
            double eoutlet = 0.5 * rhooutlet * uoutlet * uoutlet
                           + rhooutlet * lookup_U(T_clamped, P_clamped);

            if (poutlet > PRESSURE_R) {
                rho_flux[Nfaces-1] = rhooutlet * uoutlet;
                e_flux[Nfaces-1]   = uoutlet * (eoutlet + poutlet);
                mom_flux[Nfaces-1] = rhooutlet * uoutlet * uoutlet + poutlet;
            }
        } else if (u_in[N-2] == 0.0) {
            double uoutlet = 2.0 / (1.0 + GAMMA) * a_ideal;
            double poutlet = p_in[N-2] * pow(2.0 / (1.0 + GAMMA), 2.0 * GAMMA / (GAMMA - 1.0));
            double rhooutlet = GAMMA * poutlet / (uoutlet * uoutlet);
            double Toutlet = poutlet / (rhooutlet * RGAS);
            double T_clamped = fmax(Toutlet, 150.0);
            double P_clamped = fmax(poutlet, 1e5);
            double eoutlet = 0.5 * rhooutlet * uoutlet * uoutlet
                           + rhooutlet * lookup_U(T_clamped, P_clamped);

            rho_flux[Nfaces-1] = rhooutlet * uoutlet;
            e_flux[Nfaces-1]   = uoutlet * (eoutlet + poutlet);
            mom_flux[Nfaces-1] = rhooutlet * uoutlet * uoutlet + poutlet;
        }
    }

    /* Conservative update */
    /* Cell 0 (left wall BC) */
    rho_out[0] = rho_in[0] - dt / dx[0] * rho_flux[0];
    double mom1_val = rho_in[0] * u_in[0] - dt / dx[0] * (mom_flux[0] - p_in[0]);
    u_out[0] = mom1_val / rho_out[0];
    e_out[0] = e_in[0] - dt / dx[0] * e_flux[0];

    /* Interior cells */
    for (int i = 1; i < N - 1; i++) {
        rho_out[i] = rho_in[i] + dt / dx[i] * (rho_flux[i-1] - rho_flux[i]);
        double new_mom = mom[i] + dt / dx[i] * (mom_flux[i-1] - mom_flux[i]);
        u_out[i] = new_mom / rho_out[i];
        e_out[i] = e_in[i] + dt / dx[i] * (e_flux[i-1] - e_flux[i]);
    }

    /* Friction source term (implicit) */
    if (!g_no_source) {
        for (int i = 0; i < N; i++) {
            double denom_f = 1.0 + sqrt(1.0 + 2.0 * dt * fabs(u_out[i]) * lambda_fric / DIA);
            u_out[i] = 2.0 * u_out[i] / denom_f;
        }
    }

    /* Heat transfer source term */
    if (!g_no_source) {
        for (int i = 0; i < N; i++) {
            double tdiff = T0 - T_in[i];
            double deltaT = tdiff * (1.0 - exp(-dt * WATT_PER_K_PER_M3 / (rho_in[i] * Cv_in[i])));
            double deltaE = deltaT * rho_in[i] * Cv_in[i];
            e_out[i] += deltaE;
        }
    }

    /* Right boundary */
    double cfac_bc = compression_factor(T0, PRESSURE_R);
    u_out[N-1] = u_out[N-2];
    rho_out[N-1] = PRESSURE_R / (cfac_bc * RGAS * T0);
    e_out[N-1] = rho_out[N-1] * lookup_U(T0, PRESSURE_R)
               + 0.5 * rho_out[N-1] * u_out[N-1] * u_out[N-1];
}

/* ========== Recover thermodynamic state (Newton for T, then p, a) ========== */

static void recover_state(double *rho_in, double *u_in, double *e_in,
                          double *T_out, double *p_out, double *a_out,
                          double *Cv_out, int N)
{
    for (int i = 0; i < N; i++) {
        double e_int = (e_in[i] - 0.5 * rho_in[i] * u_in[i] * u_in[i]) / rho_in[i];

        if (g_ideal_gas) {
            /* Direct: no tables, no iteration */
            double Ti = e_int / CV;
            T_out[i] = Ti;
            Cv_out[i] = CV;
            p_out[i] = rho_in[i] * RGAS * Ti;
            a_out[i] = sqrt(GAMMA * p_out[i] / rho_in[i]);
        } else {
            double p_tab = p_out[i];
            if (p_tab < 1e5) p_tab = 1e5;

            double Ti = T_out[i];
            for (int iter = 0; iter < 3; iter++) {
                double cv = lookup_Cv(Ti, p_tab);
                Cv_out[i] = cv;
                Ti += (e_int - lookup_U(Ti, p_tab)) / cv;
            }
            if (Ti < 150.0) Ti = 150.0;
            T_out[i] = Ti;

            double cfac = compression_factor(Ti, p_out[i]);
            p_out[i] = Ti * rho_in[i] * RGAS * cfac;

            a_out[i] = thermo_soundspeed(Ti, p_out[i], rho_in[i], Cv_out[i]);
        }
    }
}

/* ========== Main ========== */

int main(int argc, char **argv)
{
    int use_lf = 0;
    int no_source = 0;
    int ideal_gas = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--lf") == 0) use_lf = 1;
        if (strcmp(argv[i], "--nosrc") == 0) no_source = 1;
        if (strcmp(argv[i], "--ideal") == 0) ideal_gas = 1;
    }
    g_no_source = no_source;
    g_ideal_gas = ideal_gas;

    double lambda_fric = compute_lambda();

    /* ===== Grid setup (same as MATLAB) ===== */
    int n_faces = 1;
    x_face[0] = 0.0;
    double dx_c = -DX_INIT;

    while (x_face[n_faces - 1] > -PIPE_LENGTH) {
        x_face[n_faces] = x_face[n_faces - 1] + dx_c;
        n_faces++;
        if (n_faces > MAX_CELLS) { fprintf(stderr, "Too many cells\n"); return 1; }
        if (dx_c > -DX_MAX) dx_c *= FAC;
    }
    x_face[n_faces - 1] = -PIPE_LENGTH;

    /* Reverse array so x_face goes from -pipe_length to 0 */
    for (int i = 0; i < n_faces / 2; i++) {
        double tmp = x_face[i];
        x_face[i] = x_face[n_faces - 1 - i];
        x_face[n_faces - 1 - i] = tmp;
    }

    int N = n_faces - 1;  /* number of cells */
    printf("Mode: %s%s%s\n",
           use_lf ? "Lax-Friedrich (forward Euler)" : "MUSCL-HLLC (SSP-RK2)",
           ideal_gas ? " [IDEAL GAS]" : "",
           no_source ? " [NO SOURCE]" : "");
    printf("Grid: %d cells, %d faces\n", N, n_faces);

    for (int i = 0; i < N; i++) {
        xc[i] = 0.5 * (x_face[i] + x_face[i+1]);
        dx[i] = x_face[i+1] - x_face[i];
    }

    /* ===== Initial conditions ===== */
    double cfacR = compression_factor(T0, PRESSURE_R);
    double cfacL = compression_factor(T0, PRESSURE_L);
    double densityR = PRESSURE_R / (T0 * RGAS * cfacR);
    double densityL = PRESSURE_L / (T0 * RGAS * cfacL);

    for (int i = 0; i < N; i++) {
        if (i < N - 1) {
            rho[i] = densityL;
            p[i]   = PRESSURE_L;
        } else {
            rho[i] = densityR;
            p[i]   = PRESSURE_R;
        }
        u[i] = 0.0;
        T_arr[i] = T0;
        e[i] = rho[i] * lookup_U(T0, p[i]);
    }

    /* ===== Output setup ===== */
    double dt_tseries = END_TIME / N_TSERIES;
    double dt_profiles = END_TIME / N_PROFILES;
    double next_tseries = 0.0;
    double next_profiles = 0.0;

    const char *ts_file = use_lf ? "timeseries_lf.dat" : "timeseries.dat";
    const char *pr_file = use_lf ? "profiles_lf.dat"   : "profiles.dat";
    FILE *fp_ts = fopen(ts_file, "w");
    FILE *fp_pr = fopen(pr_file, "w");
    if (!fp_ts || !fp_pr) { fprintf(stderr, "Cannot open output files\n"); return 1; }
    fprintf(fp_ts, "# t massflow u_out rho_out p_out T_out totmass p_base p_100\n");

    /* ===== Time integration ===== */
    double sim_time = 0.0;
    int nstep = 0;
    struct timespec ts_start, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    /* Initialise a_arr for first CFL */
    for (int i = 0; i < N; i++) {
        a_arr[i] = thermo_soundspeed(T_arr[i], p[i], rho[i], CV);
    }

    while (sim_time <= END_TIME) {
        nstep++;

        /* Recover thermodynamic state from conserved variables */
        recover_state(rho, u, e, T_arr, p, a_arr, Cv_real, N);

        /* CFL timestep */
        double dt = 1e30;
        for (int i = 0; i < N; i++) {
            double local_dt = CFL_NUM / (a_arr[i] + fabs(u[i])) * dx[i];
            if (local_dt < dt) dt = local_dt;
        }
        sim_time += dt;

        if (use_lf) {
            /* ===== Lax-Friedrich: single forward Euler step ===== */
            compute_lf_step(rho, u, e, p, a_arr, T_arr, Cv_real,
                            N, dt, lambda_fric,
                            rho1, u1, e1);
            memcpy(rho, rho1, N * sizeof(double));
            memcpy(u,   u1,   N * sizeof(double));
            memcpy(e,   e1,   N * sizeof(double));
        } else {
            /* ===== HLLC with SSP-RK2 ===== */

            /* Save stage-0 state */
            memcpy(rho_0, rho, N * sizeof(double));
            for (int i = 0; i < N; i++) mom_0[i] = rho[i] * u[i];
            memcpy(e_0, e, N * sizeof(double));

            /* Stage 1 */
            compute_hllc_flux(rho, u, e, p, a_arr, T_arr, Cv_real,
                              N, dt, lambda_fric, 0,
                              rho1, u1, e1);

            /* Recover state after stage 1 */
            memcpy(T1, T_arr, N * sizeof(double));
            memcpy(p1_arr, p, N * sizeof(double));
            recover_state(rho1, u1, e1, T1, p1_arr, a1_arr, Cv_real1, N);

            /* Stage 2 */
            compute_hllc_flux(rho1, u1, e1, p1_arr, a1_arr, T1, Cv_real1,
                              N, dt, lambda_fric, 1,
                              rho2, u2_arr, e2);

            /* SSP-RK2 combination: 0.5*stage0 + 0.5*stage2 */
            for (int i = 0; i < N; i++) {
                rho[i] = 0.5 * rho_0[i] + 0.5 * rho2[i];
                e[i]   = 0.5 * e_0[i]   + 0.5 * e2[i];
                double new_mom = 0.5 * mom_0[i] + 0.5 * (rho2[i] * u2_arr[i]);
                u[i] = new_mom / rho[i];
            }

            /* Re-apply BCs after averaging */
            double cfac_bc = compression_factor(T0, PRESSURE_R);
            u[N-1] = u[N-2];
            rho[N-1] = PRESSURE_R / (cfac_bc * RGAS * T0);
            e[N-1] = rho[N-1] * lookup_U(T0, PRESSURE_R)
                    + 0.5 * rho[N-1] * u[N-1] * u[N-1];
        }

        /* ===== Output ===== */
        if (sim_time >= next_tseries) {
            next_tseries += dt_tseries;
            double totmass = 0.0;
            for (int i = 0; i < N; i++) totmass += rho[i] * dx[i] * AREA;
            double mf = u[N-2] * rho[N-2];
            fprintf(fp_ts, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
                    sim_time, mf, u[N-2], rho[N-2], p[N-2], T_arr[N-2], totmass, p[0], p[100]);
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
            printf("t = %.4f / %.1f  steps = %d  steps/s = %.1f\n",
                   sim_time, END_TIME, nstep, nstep / elapsed);
        }
    }

    fclose(fp_ts);
    fclose(fp_pr);

    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    double elapsed = (ts_now.tv_sec - ts_start.tv_sec)
                   + (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
    printf("\nDone: %d steps in %.2f seconds (%.0f steps/s)\n", nstep, elapsed, nstep / elapsed);

    /* Mass conservation check */
    double totmass = 0.0;
    for (int i = 0; i < N; i++) totmass += rho[i] * dx[i] * AREA;
    printf("Final total mass: %.6f kg\n", totmass);

    return 0;
}
