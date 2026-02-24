%% 1D Shock Tube using HLLC with MUSCL-MC reconstruction and SSP-RK2
%close all
clearvars -except pressure_plot_save u_plot_save T_plot_save e_plot_save rho_plot_save
clc
%% Initialization of Parameters

Cv = 1647;
Cp = 2160;
Rgas = 8314.46 / 16.38;  % specific gas constant, MW = 16.38 kg/kmol
gamma = Cp / Cv;



dia = 1.153; %pipe inner diameter
rough = 0.017 * 1e-3; %NS empirical with heat transfer
%rough = 0.045 * 1e-3; %NS empirical without heat transfer
%rough = 0.0015 * 1e-3; %plastic
lambda = 1 / (-2*log10(rough/(3.7 * dia)))^2; %darcy friction factor (-)

area = dia.^2/4*pi;
circum = dia * pi;
Kconcrete = 4; %thermal conductivity (W/(mK))
Ksteel = 50;

thick=4.3*2.54*0.01; %thickness of concrete wall
wattperkelvinpermeter = Kconcrete/thick*circum; % W/(mK)
wattperkelvinperm3 = wattperkelvinpermeter / area; %(W/(m^3K))


% SWITCHES ================================
save_ = 1;
load_ = 0;
shocktube = 1;
model_temperature = 1; %model heat conduction from pipe walls (0 = adiabatic)
use_fric = 1; %include friction
%===========================

% Discretization and iteration time parameters ===============
  
endTime = 20;             % Physical End Time in seconds
CFL     = 0.6;
dx_init = 1;
                   %length of control volume at opening
dx_max = 0.5; %maximum length of control volume
fac = 1; %increase length of control volumes by this factor away from opening

%  Pipeline parameters   =============================
  
ns1_len = 1224 * 1000;
ns2_len = 1230 * 1000;
ns2d = 153.6 * 1000;  % to pig-trap
ns1d = 217.7 * 1000;
ns1s = 224 * 1000;
pipe_length = 4400;

pressureR   = 7.1 * 100000;  %water pressure
pressureL   = 165 * 100000; %ns1 pressure, ns2 was 103, ns1 was 165
T0 = 282;  %initial temperature of gas inside pipeline
   
%==========================================================
  
dtplot = endTime / 1000;
dt_pplot = endTime / 100;
dx_c = -dx_init;
pos = 0;
x = zeros(1);

while (x(end) > -pipe_length)

    x = [x, x(end) + dx_c];
    if (dx_c > -dx_max)
        dx_c = dx_c * fac;
    end
end

x(end) = -pipe_length;
x = x(end:-1:1);
N = length(x) - 1;
xc = 0.5 * (x(1:end-1) + x(2:end));
dx = diff(x);
time = 0;

cfacR = methane_compression_factor(T0, pressureR);
cfacL = methane_compression_factor(T0, pressureL);
    
densityR    = pressureR / (T0 * Rgas) / cfacR;       densityL    = pressureL / (T0 * Rgas) / cfacL;
    
cl = sqrt(gamma*pressureL/densityL);
umid = 2 / (1 + gamma) * cl;
fact = 1. - 0.5 * (gamma-1) * umid / cl;
    
rhomid = densityL * fact ^ (2. / (gamma - 1));
    
pmid = pressureL * fact ^ (2. * gamma / (gamma-1));
    
rho     = zeros(N,1) ;
p       = zeros(N,1) ;
u       = zeros(N,1) ;
    
%PV=nRT => PV ~ rho*T
    
    for i =  1:N
        if (i < N)
            rho(i)  = densityL  ;
            p(i)    = pressureL ;
        else
            rho(i)  = densityR  ;
            p(i)    = pressureR ;
        end
    end
%%==================initialize a pulse =============================
	%plength = 120; %initial length of pulse
	n_cells = 9; %floor(plength/dx_init); 
        %mid=floor(N/2);
        %mid=N-300-n_cells;
        %rho(mid:1:mid+n_cells) = densityL*6;  %adiabatic compression of 1m3 methane: V=1/3, P=P0x5.43 and T=T0x1.404
        %p(mid:1:mid+n_cells)=pressureL*14;
%%==================end initialize a pulse =============================
            
    if (load_)
        load('workspace.mat');
        time = 0;
        inc  = ceil(200 / dx);
        rho(inc:N/2) = rho(1:N/2 - inc + 1);
        rho(1:inc - 1) = densityL;
        rho(N/2+1:end) = densityR;
        rho(N/2+1:end) = densityR;

    
        p(inc:N/2) = p(1:N/2 - inc + 1);
        p(1:inc - 1) = pressureL;
        p(N/2+1:end) = pressureR;
        u(inc:N/2) = u(1:N/2 - inc + 1);
        u(1:inc - 1) = 0;
        u(N/2+1:end) = 0;
    
    end

u0 = u;
p0 = p;
rho0 = rho;


T = ones(N,1) * T0; % Erik: T0 should be a function of pipe length. It was average 279 for northern end but 285 at outlet

% temperature of pulse =========================
		%	T(mid:1:mid+n_cells)= 479; %T0*1.424; %result of adiabatic compression to third of volume. 13.6MJ
%======================================

e = rho .* methane_U(T, p);  % real gas internal energy from CoolProp (u=0 so no KE)


%%
new_rho = rho ;
new_u   = u   ;
new_e   = e   ;
new_T = T;
new_mom = rho.*u;

totmass = [];
mass_right = [];
mass_left = [];
e_right = [];
e_left = [];
plot_times = [];
tot_e = [];
p_base = [];
massflow_out = [];
u_out = [];
rho_out = [];
e_out = [];
T_out = [];
p_out = [];

pressure_plot = [];
e_plot = [];
u_plot = [];
rho_plot = [];
T_plot = [];
pplot_times = [];

plot_time = time;

pplot_time = time;

nstep = 0;

tstart = tic;

while time <= endTime
    

    nstep = nstep + 1;
    
    % Recover T from conserved energy using real gas internal energy
    % Newton iteration: solve methane_U(T, p) = e_int_per_mass for T
    u_target = (e - 0.5*rho.*u.^2) ./ rho;  % target specific internal energy
    p_tab = max(p, 1e5);  % clamp pressure for table lookups (table min = 1 bar)
    for newton_iter = 1:3
        Cv_real = methane_Cv(T, p_tab);
        T = T + (u_target - methane_U(T, p_tab)) ./ Cv_real;
    end
    T = max(T, 150);

    cfac = methane_compression_factor(T, p);
    p = T .* rho * Rgas .* cfac;

    p_tab = max(p, 1e5);
    a = methane_soundspeed(T, p_tab);  % real gas sound speed from CoolProp
        
    dt      = min(CFL./(a(:)+abs(u(:))) .* dx(:)) ;  % adjustable Time Step
    time    = time + dt ;

    %fprintf("time: %f\n", time);
    % ============= SSP-RK2 Stage 1: Forward Euler predictor =============
    mom = rho .* u;

    % --- MUSCL reconstruction with MC limiter ---
    slope_rho = zeros(N,1);
    slope_mom = zeros(N,1);
    slope_e   = zeros(N,1);
    for ii = 2:N-1
        dx_L = xc(ii) - xc(ii-1);
        dx_R = xc(ii+1) - xc(ii);
        slope_rho(ii) = mc_limiter((rho(ii)-rho(ii-1))/dx_L, (rho(ii+1)-rho(ii))/dx_R);
        slope_mom(ii) = mc_limiter((mom(ii)-mom(ii-1))/dx_L, (mom(ii+1)-mom(ii))/dx_R);
        slope_e(ii)   = mc_limiter((e(ii)-e(ii-1))/dx_L,     (e(ii+1)-e(ii))/dx_R);
    end

    % Reconstruct left/right states at each face i+1/2 (face between cell i and i+1)
    rhoL_f = rho(1:end-1) + 0.5*dx(1:end-1)'.*slope_rho(1:end-1);
    rhoR_f = rho(2:end)   - 0.5*dx(2:end)'  .*slope_rho(2:end);
    momL_f = mom(1:end-1) + 0.5*dx(1:end-1)'.*slope_mom(1:end-1);
    momR_f = mom(2:end)   - 0.5*dx(2:end)'  .*slope_mom(2:end);
    eL_f   = e(1:end-1)   + 0.5*dx(1:end-1)'.*slope_e(1:end-1);
    eR_f   = e(2:end)     - 0.5*dx(2:end)'  .*slope_e(2:end);

    % Enforce positivity on reconstructed density and energy
    rhoL_f = max(rhoL_f, 1e-6);
    rhoR_f = max(rhoR_f, 1e-6);
    eL_f   = max(eL_f, 1e-6);
    eR_f   = max(eR_f, 1e-6);

    uL_f = momL_f ./ rhoL_f;
    uR_f = momR_f ./ rhoR_f;
    % Pressure from reconstructed states: use cell-center pressure as approximation
    pL_f = p(1:end-1);
    pR_f = p(2:end);
    % Sound speed: use cell-center values
    aL_f = a(1:end-1);
    aR_f = a(2:end);

    % --- HLLC flux ---
    SL = min(uL_f - aL_f, uR_f - aR_f);
    SR = max(uL_f + aL_f, uR_f + aR_f);
    denom = rhoL_f.*(SL - uL_f) - rhoR_f.*(SR - uR_f);
    denom(abs(denom) < 1e-30) = 1e-30;  % avoid division by zero
    SM = (pR_f - pL_f + rhoL_f.*uL_f.*(SL - uL_f) - rhoR_f.*uR_f.*(SR - uR_f)) ./ denom;

    % Physical fluxes from left and right states
    FL_rho = momL_f;
    FL_mom = rhoL_f.*uL_f.*uL_f + pL_f;
    FL_e   = uL_f .* (eL_f + pL_f);
    FR_rho = momR_f;
    FR_mom = rhoR_f.*uR_f.*uR_f + pR_f;
    FR_e   = uR_f .* (eR_f + pR_f);

    % Star states (left)
    dSL = SL - uL_f;
    dSL(abs(dSL) < 1e-30) = 1e-30;
    rho_starL = rhoL_f .* (SL - uL_f) ./ (SL - SM);
    mom_starL = rho_starL .* SM;
    E_starL   = rho_starL .* (eL_f./rhoL_f + (SM - uL_f).*(SM + pL_f./(rhoL_f.*dSL)));

    % Star states (right)
    dSR = SR - uR_f;
    dSR(abs(dSR) < 1e-30) = 1e-30;
    rho_starR = rhoR_f .* (SR - uR_f) ./ (SR - SM);
    mom_starR = rho_starR .* SM;
    E_starR   = rho_starR .* (eR_f./rhoR_f + (SM - uR_f).*(SM + pR_f./(rhoR_f.*dSR)));

    % Assemble HLLC flux based on wave pattern
    rho_flux_e = zeros(N-1, 1);
    mom_flux_e = zeros(N-1, 1);
    energy_flux_e = zeros(N-1, 1);

    idx1 = (SL >= 0);
    rho_flux_e(idx1) = FL_rho(idx1);
    mom_flux_e(idx1) = FL_mom(idx1);
    energy_flux_e(idx1) = FL_e(idx1);

    idx2 = (~idx1) & (SM >= 0);
    rho_flux_e(idx2) = FL_rho(idx2) + SL(idx2).*(rho_starL(idx2) - rhoL_f(idx2));
    mom_flux_e(idx2) = FL_mom(idx2) + SL(idx2).*(mom_starL(idx2) - momL_f(idx2));
    energy_flux_e(idx2) = FL_e(idx2) + SL(idx2).*(E_starL(idx2) - eL_f(idx2));

    idx3 = (~idx1) & (~idx2) & (SR >= 0);
    rho_flux_e(idx3) = FR_rho(idx3) + SR(idx3).*(rho_starR(idx3) - rhoR_f(idx3));
    mom_flux_e(idx3) = FR_mom(idx3) + SR(idx3).*(mom_starR(idx3) - momR_f(idx3));
    energy_flux_e(idx3) = FR_e(idx3) + SR(idx3).*(E_starR(idx3) - eR_f(idx3));

    idx4 = (~idx1) & (~idx2) & (~idx3);
    rho_flux_e(idx4) = FR_rho(idx4);
    mom_flux_e(idx4) = FR_mom(idx4);
    energy_flux_e(idx4) = FR_e(idx4);

    choked_bc = 1;

    if (choked_bc)
        % Use ideal gas sound speed for choked BC (consistent with
        % ideal gas isentropic relations used here)
        a_ideal = sqrt(gamma * p(N-1) / rho(N-1));

        if (u(N-1) < a_ideal && u(N - 1) > 0)

            vc = u(N - 1) / a_ideal;

            xd = vc * (1 + gamma) / (2 + vc * (gamma - 1));

            uoutlet = u(N - 1) / xd;

            c0 = 0.5 * (1 + gamma) * uoutlet;

            p0 = p(N - 1) * (c0 / a_ideal) ^ (2*gamma/ (gamma - 1));

            poutlet = p0 * (2 / (1 + gamma))^(2*gamma/(gamma - 1));
            rhooutlet = gamma * poutlet / uoutlet^2;
            Toutlet = poutlet / (rhooutlet * Rgas);  % ideal gas T (Z=1 in this context)
            eoutlet = uoutlet.^2*rhooutlet/2 + rhooutlet * methane_U(max(Toutlet,150), max(poutlet, 1e5));

            if (poutlet > pressureR)
                rho_flux_e(end) = rhooutlet * uoutlet;
                energy_flux_e(end) = uoutlet * (eoutlet + poutlet);
                mom_flux_e(end) = rhooutlet * uoutlet^2 + poutlet;
            end

        elseif (u(N - 1) == 0)
            uoutlet = 2 / (1 + gamma) * a_ideal;
            poutlet = p(N - 1)*(2 / (1 + gamma)) ^ (2 * gamma / (gamma - 1));
            rhooutlet = gamma * poutlet / uoutlet^2;

            Toutlet = poutlet / (rhooutlet * Rgas);  % ideal gas T (Z=1)
            eoutlet = uoutlet.^2*rhooutlet/2 + rhooutlet * methane_U(max(Toutlet,150), max(poutlet, 1e5));

            rho_flux_e(end) = rhooutlet * uoutlet;
            energy_flux_e(end) = uoutlet * (eoutlet + poutlet);
            mom_flux_e(end) = rhooutlet * uoutlet^2 + poutlet;
        end
    end

    % Stage 1 conservative update
    new_rho(2:end-1) = rho(2:end-1) + dt ./ dx(2:end-1)' .* (rho_flux_e(1:end-1) - rho_flux_e(2:end));
    new_mom(2:end-1) = mom(2:end-1) + dt ./ dx(2:end-1)' .* (mom_flux_e(1:end-1) - mom_flux_e(2:end));
    new_e(2:end-1) = e(2:end-1) + dt./dx(2:end-1)' .* (energy_flux_e(1:end-1) - energy_flux_e(2:end));

    new_rho(1) = rho(1) - dt/dx(1) * rho_flux_e(1);
    mom1_stage = rho(1)*u(1) - dt/dx(1) * (mom_flux_e(1) - (0*rho(1)*u(1).^2 +p(1)));
    new_u(1) = mom1_stage / new_rho(1);
    new_e(1) = e(1) - dt/dx(1) * (energy_flux_e(1) - 0*(u(1) * ( e(1) + p(1) )));

    new_u(2:end-1) = new_mom(2:end-1) ./ new_rho(2:end-1);

    if (use_fric)
        u_fric = (2 * new_u) ./ (1+sqrt(1+(2*dt.*abs(new_u) * lambda) / dia));
        new_u = u_fric;
    end

    if (model_temperature)
        tdiff = T0 - T;
        deltaT = tdiff .* (1 - exp(-dt * wattperkelvinperm3 ./ (rho .* Cv_real)));
        deltaE = deltaT .* rho .* Cv_real;
        new_e = new_e + deltaE;
    end

    new_u(N) = new_u(N - 1);
    cfac_bc = methane_compression_factor(T0, pressureR);
    new_rho(N) = pressureR / (cfac_bc * Rgas * T0);
    new_e(N) = new_rho(N) * methane_U(T0, pressureR) + 0.5*new_rho(N) .* new_u(N) .* new_u(N);

    % Save stage-0 state for SSP-RK2 averaging
    rho_0 = rho;
    mom_0 = mom;
    e_0   = e;

    % Promote stage 1 to current state
    rho1 = new_rho;
    u1   = new_u;
    e1   = new_e;
    mom1 = rho1 .* u1;
    T1   = T;  % will be updated below

    % Recover thermodynamic state for stage 1
    u_target1 = (e1 - 0.5*rho1.*u1.^2) ./ rho1;
    p_tab1 = max(p, 1e5);
    for newton_iter = 1:3
        Cv_real1 = methane_Cv(T1, p_tab1);
        T1 = T1 + (u_target1 - methane_U(T1, p_tab1)) ./ Cv_real1;
    end
    T1 = max(T1, 150);
    cfac1 = methane_compression_factor(T1, p);
    p1 = T1 .* rho1 * Rgas .* cfac1;
    p_tab1 = max(p1, 1e5);
    a1 = methane_soundspeed(T1, p_tab1);

    % ============= SSP-RK2 Stage 2: second Euler step from stage 1 =============
    % --- MUSCL reconstruction on stage-1 state ---
    slope_rho2 = zeros(N,1);
    slope_mom2 = zeros(N,1);
    slope_e2   = zeros(N,1);
    for ii = 2:N-1
        dx_L = xc(ii) - xc(ii-1);
        dx_R = xc(ii+1) - xc(ii);
        slope_rho2(ii) = mc_limiter((rho1(ii)-rho1(ii-1))/dx_L, (rho1(ii+1)-rho1(ii))/dx_R);
        slope_mom2(ii) = mc_limiter((mom1(ii)-mom1(ii-1))/dx_L, (mom1(ii+1)-mom1(ii))/dx_R);
        slope_e2(ii)   = mc_limiter((e1(ii)-e1(ii-1))/dx_L,     (e1(ii+1)-e1(ii))/dx_R);
    end

    rhoL_f2 = rho1(1:end-1) + 0.5*dx(1:end-1)'.*slope_rho2(1:end-1);
    rhoR_f2 = rho1(2:end)   - 0.5*dx(2:end)'  .*slope_rho2(2:end);
    momL_f2 = mom1(1:end-1) + 0.5*dx(1:end-1)'.*slope_mom2(1:end-1);
    momR_f2 = mom1(2:end)   - 0.5*dx(2:end)'  .*slope_mom2(2:end);
    eL_f2   = e1(1:end-1)   + 0.5*dx(1:end-1)'.*slope_e2(1:end-1);
    eR_f2   = e1(2:end)     - 0.5*dx(2:end)'  .*slope_e2(2:end);

    rhoL_f2 = max(rhoL_f2, 1e-6);
    rhoR_f2 = max(rhoR_f2, 1e-6);
    eL_f2   = max(eL_f2, 1e-6);
    eR_f2   = max(eR_f2, 1e-6);

    uL_f2 = momL_f2 ./ rhoL_f2;
    uR_f2 = momR_f2 ./ rhoR_f2;
    pL_f2 = p1(1:end-1);
    pR_f2 = p1(2:end);
    aL_f2 = a1(1:end-1);
    aR_f2 = a1(2:end);

    % --- HLLC flux (stage 2) ---
    SL2 = min(uL_f2 - aL_f2, uR_f2 - aR_f2);
    SR2 = max(uL_f2 + aL_f2, uR_f2 + aR_f2);
    denom2 = rhoL_f2.*(SL2 - uL_f2) - rhoR_f2.*(SR2 - uR_f2);
    denom2(abs(denom2) < 1e-30) = 1e-30;
    SM2 = (pR_f2 - pL_f2 + rhoL_f2.*uL_f2.*(SL2 - uL_f2) - rhoR_f2.*uR_f2.*(SR2 - uR_f2)) ./ denom2;

    FL_rho2 = momL_f2;
    FL_mom2 = rhoL_f2.*uL_f2.*uL_f2 + pL_f2;
    FL_e2   = uL_f2 .* (eL_f2 + pL_f2);
    FR_rho2 = momR_f2;
    FR_mom2 = rhoR_f2.*uR_f2.*uR_f2 + pR_f2;
    FR_e2   = uR_f2 .* (eR_f2 + pR_f2);

    dSL2 = SL2 - uL_f2;
    dSL2(abs(dSL2) < 1e-30) = 1e-30;
    rho_sL2 = rhoL_f2 .* (SL2 - uL_f2) ./ (SL2 - SM2);
    mom_sL2 = rho_sL2 .* SM2;
    E_sL2   = rho_sL2 .* (eL_f2./rhoL_f2 + (SM2 - uL_f2).*(SM2 + pL_f2./(rhoL_f2.*dSL2)));

    dSR2 = SR2 - uR_f2;
    dSR2(abs(dSR2) < 1e-30) = 1e-30;
    rho_sR2 = rhoR_f2 .* (SR2 - uR_f2) ./ (SR2 - SM2);
    mom_sR2 = rho_sR2 .* SM2;
    E_sR2   = rho_sR2 .* (eR_f2./rhoR_f2 + (SM2 - uR_f2).*(SM2 + pR_f2./(rhoR_f2.*dSR2)));

    rho_flux_e2 = zeros(N-1, 1);
    mom_flux_e2 = zeros(N-1, 1);
    energy_flux_e2 = zeros(N-1, 1);

    idx1 = (SL2 >= 0);
    rho_flux_e2(idx1) = FL_rho2(idx1);
    mom_flux_e2(idx1) = FL_mom2(idx1);
    energy_flux_e2(idx1) = FL_e2(idx1);

    idx2 = (~idx1) & (SM2 >= 0);
    rho_flux_e2(idx2) = FL_rho2(idx2) + SL2(idx2).*(rho_sL2(idx2) - rhoL_f2(idx2));
    mom_flux_e2(idx2) = FL_mom2(idx2) + SL2(idx2).*(mom_sL2(idx2) - momL_f2(idx2));
    energy_flux_e2(idx2) = FL_e2(idx2) + SL2(idx2).*(E_sL2(idx2) - eL_f2(idx2));

    idx3 = (~idx1) & (~idx2) & (SR2 >= 0);
    rho_flux_e2(idx3) = FR_rho2(idx3) + SR2(idx3).*(rho_sR2(idx3) - rhoR_f2(idx3));
    mom_flux_e2(idx3) = FR_mom2(idx3) + SR2(idx3).*(mom_sR2(idx3) - momR_f2(idx3));
    energy_flux_e2(idx3) = FR_e2(idx3) + SR2(idx3).*(E_sR2(idx3) - eR_f2(idx3));

    idx4 = (~idx1) & (~idx2) & (~idx3);
    rho_flux_e2(idx4) = FR_rho2(idx4);
    mom_flux_e2(idx4) = FR_mom2(idx4);
    energy_flux_e2(idx4) = FR_e2(idx4);

    % Choked BC for stage 2
    if (choked_bc)
        a_ideal2 = sqrt(gamma * p1(N-1) / rho1(N-1));

        if (u1(N-1) < a_ideal2 && u1(N - 1) > 0)
            vc2 = u1(N - 1) / a_ideal2;
            xd2 = vc2 * (1 + gamma) / (2 + vc2 * (gamma - 1));
            uoutlet2 = u1(N - 1) / xd2;
            c02 = 0.5 * (1 + gamma) * uoutlet2;
            p02 = p1(N - 1) * (c02 / a_ideal2) ^ (2*gamma/ (gamma - 1));
            poutlet2 = p02 * (2 / (1 + gamma))^(2*gamma/(gamma - 1));
            rhooutlet2 = gamma * poutlet2 / uoutlet2^2;
            Toutlet2 = poutlet2 / (rhooutlet2 * Rgas);
            eoutlet2 = uoutlet2.^2*rhooutlet2/2 + rhooutlet2 * methane_U(max(Toutlet2,150), max(poutlet2, 1e5));

            if (poutlet2 > pressureR)
                rho_flux_e2(end) = rhooutlet2 * uoutlet2;
                energy_flux_e2(end) = uoutlet2 * (eoutlet2 + poutlet2);
                mom_flux_e2(end) = rhooutlet2 * uoutlet2^2 + poutlet2;
            end

        elseif (u1(N - 1) == 0)
            uoutlet2 = 2 / (1 + gamma) * a_ideal2;
            poutlet2 = p1(N - 1)*(2 / (1 + gamma)) ^ (2 * gamma / (gamma - 1));
            rhooutlet2 = gamma * poutlet2 / uoutlet2^2;
            Toutlet2 = poutlet2 / (rhooutlet2 * Rgas);
            eoutlet2 = uoutlet2.^2*rhooutlet2/2 + rhooutlet2 * methane_U(max(Toutlet2,150), max(poutlet2, 1e5));

            rho_flux_e2(end) = rhooutlet2 * uoutlet2;
            energy_flux_e2(end) = uoutlet2 * (eoutlet2 + poutlet2);
            mom_flux_e2(end) = rhooutlet2 * uoutlet2^2 + poutlet2;
        end
    end

    % Stage 2 conservative update (from stage 1 state)
    new_rho2 = rho1;
    new_mom2 = mom1;
    new_e2   = e1;
    new_u2   = u1;

    new_rho2(2:end-1) = rho1(2:end-1) + dt ./ dx(2:end-1)' .* (rho_flux_e2(1:end-1) - rho_flux_e2(2:end));
    new_mom2(2:end-1) = mom1(2:end-1) + dt ./ dx(2:end-1)' .* (mom_flux_e2(1:end-1) - mom_flux_e2(2:end));
    new_e2(2:end-1)   = e1(2:end-1)   + dt./dx(2:end-1)' .* (energy_flux_e2(1:end-1) - energy_flux_e2(2:end));

    new_rho2(1) = rho1(1) - dt/dx(1) * rho_flux_e2(1);
    mom1_stage2 = rho1(1)*u1(1) - dt/dx(1) * (mom_flux_e2(1) - (0*rho1(1)*u1(1).^2 +p1(1)));
    new_u2(1) = mom1_stage2 / new_rho2(1);
    new_e2(1) = e1(1) - dt/dx(1) * (energy_flux_e2(1) - 0*(u1(1) * ( e1(1) + p1(1) )));

    new_u2(2:end-1) = new_mom2(2:end-1) ./ new_rho2(2:end-1);

    if (use_fric)
        u_fric2 = (2 * new_u2) ./ (1+sqrt(1+(2*dt.*abs(new_u2) * lambda) / dia));
        new_u2 = u_fric2;
    end

    if (model_temperature)
        tdiff1 = T0 - T1;
        deltaT1 = tdiff1 .* (1 - exp(-dt * wattperkelvinperm3 ./ (rho1 .* Cv_real1)));
        deltaE1 = deltaT1 .* rho1 .* Cv_real1;
        new_e2 = new_e2 + deltaE1;
    end

    new_u2(N) = new_u2(N - 1);
    new_rho2(N) = pressureR / (cfac_bc * Rgas * T0);
    new_e2(N) = new_rho2(N) * methane_U(T0, pressureR) + 0.5*new_rho2(N) .* new_u2(N) .* new_u2(N);

    % ============= SSP-RK2 combination: 0.5*stage0 + 0.5*stage2 =============
    new_rho = 0.5*rho_0 + 0.5*new_rho2;
    new_e   = 0.5*e_0   + 0.5*new_e2;
    new_mom = 0.5*mom_0  + 0.5*(new_rho2 .* new_u2);
    new_u   = new_mom ./ new_rho;

    % Re-apply BCs after averaging
    new_u(N) = new_u(N - 1);
    new_rho(N) = pressureR / (cfac_bc * Rgas * T0);
    new_e(N) = new_rho(N) * methane_U(T0, pressureR) + 0.5*new_rho(N) .* new_u(N) .* new_u(N);
    
    if (time >= plot_time)
        fprintf("time: %f\n", time);
        telapsed = toc(tstart);
        fprintf("steps per second: %f\n", nstep / telapsed);
        plot_time = plot_time + dtplot;
        plot_times = [plot_times time];
        totmass = [totmass sum(rho(:).*dx(:).*area)];
        mass_right = [mass_right sum(rho(floor(N/2) + 1:end).*dx(floor(N/2) + 1:end)'.*area)];
        mass_left = [mass_left sum(rho(1:floor(N/2)-1).*dx(1:floor(N/2)-1)'.*area)];
        tot_e = [tot_e sum(e)];
        e_left = [e_left sum(e(1:floor(N/2)-1))];
        e_right = [e_right sum(e(floor(N/2):end))];
        p_base = [p_base p(1)];
        %massflow_out = [massflow_out rho_flux_e(end)];
        massflow_out = [massflow_out u(N - 1) * rho(N - 1)];
        u_out = [u_out u(N - 1)];
        rho_out = [rho_out, rho(N - 1)];
        e_out = [e_out e(N - 1)];
        T_out = [T_out T(N - 1)];
        p_out = [p_out p(N - 1)];
    end

    if (time >= pplot_time) 
        pplot_time = pplot_time + dt_pplot;
        pplot_times = [pplot_times time];
        pressure_plot = [pressure_plot p];
        e_plot = [e_plot e];
        u_plot = [u_plot u];
        rho_plot = [rho_plot rho];
        T_plot = [T_plot T];
    end
    
     rho     = new_rho ;
     u       = new_u ;
     e       = new_e ;
end

if (save_)
    save("workspace.mat");
    rho_plot_save = rho_plot;
    e_plot_save = e_plot;
    T_plot_save = T_plot;
    pressure_plot_save = pressure_plot;
    u_plot_save = u_plot;
end

% pressure = dlmread('pressure.dat') ;
% density  = dlmread('density.dat')  ;
% velocity = dlmread('velocity.dat') ;
% 
% figure(1)
% hold on
% %plot(density(:,1),density(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
% plot(xc, rho,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
% xlabel(' x ','FontSize',18,'FontWeight','bold');
% ylabel(' Density ','FontSize',18,'FontWeight','bold');
% legend('Lax Friedrich','Location','northeast','FontSize',16,'FontWeight','bold');
% %print(gcf,'Density.jpg','-dpng','-r300');
% 
% figure(2)
% hold on
% %plot(pressure(:,1),pressure(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
% plot(xc, p,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
% xlabel(' x ','FontSize',18,'FontWeight','bold');
% ylabel(' Pressure ','FontSize',18,'FontWeight','bold');
% legend('Lax Friedrich','Location','northeast','FontSize',16,'FontWeight','bold');
% %print(gcf,'Pressure.jpg','-dpng','-r300');
% 
% figure(3)
% hold on
% %plot(velocity(:,1),velocity(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
% plot(xc, u,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
% xlabel(' x ','FontSize',18,'FontWeight','bold');
% ylabel(' Velocity ','FontSize',18,'FontWeight','bold');
% legend('Lax Friedrich','Location','south','FontSize',16,'FontWeight','bold');
%print(gcf,'Velocity.jpg','-dpng','-r300');
