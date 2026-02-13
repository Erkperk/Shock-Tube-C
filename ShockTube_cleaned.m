%% 1D Shock Tube using Lax-Friedrich and Adjustable Time Stepping
%close all
clearvars -except pressure_plot_save u_plot_save T_plot_save e_plot_save rho_plot_save
clc
%% Initialization of Parameters

Cv = 1622;
Cp = 2125;
Rgas = 8314.46 / 16.71;  % specific gas constant, MW = 16.71 kg/kmol
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
CFL     = 0.3;            %was 0.3
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
    if (1)
        mom = rho .* u;      %momentum = density * velocity
        lambda_t = abs(u) + a;

        lambda_e = max(lambda_t(1:end-1), lambda_t(2:end));

        mom_flux_t = rho.*u.*u + p;  
        energy_flux_t = u .* (e + p);

        rho_flux_e = 0.5 * (mom(1:end-1) + mom(2:end)) + 0.5 * lambda_e .* (rho(1:end-1) - rho(2:end));
        mom_flux_e = 0.5 * (mom_flux_t(1:end-1) + mom_flux_t(2:end)) + 0.5 * lambda_e .* (mom(1:end-1) - mom(2:end));
        energy_flux_e = 0.5 * (energy_flux_t(1:end-1) + energy_flux_t(2:end)) + 0.5 * lambda_e .* (e(1:end-1) - e(2:end));

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


        new_rho(2:end-1) = rho(2:end-1) + dt ./ dx(2:end-1)' .* (rho_flux_e(1:end-1) - rho_flux_e(2:end));
        new_mom(2:end-1) = mom(2:end-1) + dt ./ dx(2:end-1)' .* (mom_flux_e(1:end-1) - mom_flux_e(2:end));
        new_e(2:end-1) = e(2:end-1) + dt./dx(2:end-1)' .* (energy_flux_e(1:end-1) - energy_flux_e(2:end));



        new_rho(1) = rho(1) - dt/dx(1) * rho_flux_e(1);
        mom1 = rho(1)*u(1) - dt/dx(1) * (mom_flux_e(1) - (0*rho(1)*u(1).^2 +p(1)));
        new_u(1) = mom1 / new_rho(1);
        new_e(1) = e(1) - dt/dx(1) * (energy_flux_e(1) - 0*(u(1) * ( e(1) + p(1) )));

        new_u(2:end-1) = new_mom(2:end-1) ./ new_rho(2:end-1);

                
        if (use_fric)
            u_fric = (2 * new_u) ./ (1+sqrt(1+(2*dt.*abs(new_u) * lambda) / dia));
            
            new_u = u_fric;
        end
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
