%% 1D Shock Tube using Lax-Friedrich and Adjustable Time Stepping
%close all
clearvars -except pressure_plot_save u_plot_save T_plot_save e_plot_save rho_plot_save
clc
%% Initialization of Parameters

Cv = 1709.;
Cp = 2232.;
Rgas = Cp - Cv;
gamma = Cp / Cv;



dia = 1.153; %pipe inner diameter
%rough = 0.05 * 1e-3; %steel
rough = 0.017 * 1e-3; %NS empirical with heat transfer
%rough = 0.045 * 1e-3; %NS empirical without heat transfer
%rough = 0.0015 * 1e-3; %plastic
lambda = 1 / (-2*log10(rough/(3.7 * dia)))^2; %darcy friction factor (-)
%lambda = 0.1; %fanning friction coefficient

area = dia.^2/4*pi;
circum = dia * pi;
Kconcrete = 4; %thermal conductivity (W/(mK))
Ksteel = 50;

thick=4.3*2.54*0.01; %thickness of concrete wall
wattperkelvinpermeter = Kconcrete/thick*circum; % W/(mK)
wattperkelvinperm3 = wattperkelvinpermeter / area; %(W/(m^3K))

save_ = 1;
load_ = 0;

shocktube = 1;
model_temperature = 1; %model heat conduction from pipe walls (0 = adiabatic)
use_fric = 1; %include friction

endTime = 400;;             % Physical End Time


dtplot = endTime / 1000;
dt_pplot = endTime / 100;
CFL     = 0.3;   %was 0.3

dx_init = 1.2; %length of control volume at opening
dx_max = 1.2; %maximum length of control volume
fac = 1.0; %increase length of control volumes by this factor away from opening


dx_c = -dx_init;
pos = 0;

x = zeros(1);

ns1_len = 1224 * 1000;
ns2_len = 1230 * 1000;

ns2d = 150 * 1000;
ns1d = 217.7 * 1000;
ns1s = 224 * 1000;


pipe_length = ns2d*2;


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

time    = 0 ;


if (shocktube)
    
  
    pressureR   = 8.1 * 100000;       
    %pressureR   = 0.1 * 100000;       
    %pressureL   = 103 * 100000; %ns2
    %

    pressureL   = 165 * 100000; %ns1

    %%"nordstream" conditions

       
   
    
    T0 = 283;
   

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
        %if i<=N/2
        if (i < N)
            rho(i)  = densityL  ;
            p(i)    = pressureL ;
        else
            rho(i)  = densityR  ;
            p(i)    = pressureR ;
        end
    end
        mid=floor(N/2);
        rho(mid:1:mid) = densityL*3;
        p(mid:1:mid)=pressureL*10;
            
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

else     %noSt shocktub
	    

	    'hiiiiiiiiiiiiiiiiiiiiiiiiiiiIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'

    T0 = 280;
    pressureL   = 165 * 100000;
    densityL    = pressureL / (T0 * Rgas);
    %pressureS = 330 * 100000;
    %pressureL / densityL^gamma == pressureS / densityS^gamma

    densityS = 2 * densityL;
    pressureS = pressureL * (densityS / densityL) ^ gamma;

    pressureR = pressureL;
    densityR = densityL;
    
    cl = sqrt(gamma*pressureL/densityL);

    %(2 * ncell + 1) * dx = 1 => (1/dx - 1) / 2

    %ncell = ceil( (1. / dx - 1) / 2);
    ncell = 100;

    rho = ones(N, 1) * densityL;
    rho(N / 2 - ncell:N/2 + ncell) = densityS;

    p = ones(N, 1) * pressureL;
    p(N / 2 - ncell:N/2 + ncell) = pressureS;

    u = zeros(N, 1);


end

u0 = u;
p0 = p;
rho0 = rho;

%T = p ./ (rho .* Rgas);

T = ones(N,1) * T0; % Erik: T0 should be a function of pipe length. It was average 279 for northern end but 285 at outlet
T(mid)=T0*3;
e = T * Cv .* rho;  % energy


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
    
    T = (e - 0.5*rho.*u.^2) ./ (Cv .* rho);
    %previous time step pressure to avoid making this insanely
    %complicated...
    cfac = methane_compression_factor(T, p);
    
    p = T .* rho * Rgas .* cfac;

    a = sqrt(cfac .* T *Rgas * gamma);
    lamda   = max(a) ;
    max_vel = max(abs(u)) ;
    
    dt      = min(CFL./(a(:)+abs(u(:))) .* dx(:)) ;  % adjustable Time Step
    time    = time + dt ;

    %fprintf("time: %f\n", time);
    if (1)
        mom = rho .* u;      %momentum = density * velocity
        lambda_t = abs(u) + sqrt(gamma * p ./ rho);

        lambda_e = max(lambda_t(1:end-1), lambda_t(2:end));

        mom_flux_t = rho.*u.*u + p;  
        energy_flux_t = u .* (e + p);

        rho_flux_e = 0.5 * (mom(1:end-1) + mom(2:end)) + 0.5 * lambda_e .* (rho(1:end-1) - rho(2:end));
        mom_flux_e = 0.5 * (mom_flux_t(1:end-1) + mom_flux_t(2:end)) + 0.5 * lambda_e .* (mom(1:end-1) - mom(2:end));
        energy_flux_e = 0.5 * (energy_flux_t(1:end-1) + energy_flux_t(2:end)) + 0.5 * lambda_e .* (e(1:end-1) - e(2:end));

        choked_bc = 1;

        if (choked_bc)
            if (u(N-1) < a(N-1) && u(N - 1) > 0)
                
                vc = u(N - 1) / a(N - 1);

                xd = vc * (1 + gamma) / (2 + vc * (gamma - 1));

                uoutlet = u(N - 1) / xd;

                c0 = 0.5 * (1 + gamma) * uoutlet;

                p0 = p(N - 1) * (c0 / a(N - 1)) ^ (2*gamma/ (gamma - 1));

                poutlet = p0 * (2 / (1 + gamma))^(2*gamma/(gamma - 1));
                rhooutlet = gamma * poutlet / uoutlet^2;
                eoutlet = uoutlet.^2*rhooutlet/2 + poutlet / (gamma - 1) / cfac(N - 1);

                Toutlet = (eoutlet - 0.5*rhooutlet.*uoutlet.^2) ./ (Cv .* rhooutlet);

                if (poutlet > pressureR)
                    rho_flux_e(end) = rhooutlet * uoutlet;
                    energy_flux_e(end) = uoutlet * (eoutlet + poutlet);
                    mom_flux_e(end) = rhooutlet * uoutlet^2 + poutlet;
                end

            elseif (u(N - 1) == 0)
                uoutlet = 2 / (1 + gamma) * a(N - 1);
                poutlet = p(N - 1)*(2 / (1 + gamma)) ^ (2 * gamma / (gamma - 1));
                rhooutlet = gamma * poutlet / uoutlet^2;

                eoutlet = uoutlet.^2*rhooutlet/2 + poutlet / (gamma - 1) ./ cfac(N-1);

                Toutlet = (eoutlet - 0.5*rhooutlet.*uoutlet.^2) ./ (Cv .* rhooutlet);

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
    
    if (0)
    for i=2:N-1
        

        mom_R       = rho(i+1)*u(i+1) ;     rho_R = rho(i+1) ;      u_R = u(i+1) ;      p_R = p(i+1) ;
        mom_P       = rho(i)*u(i)     ; 	rho_P = rho(i)   ;      u_P = u(i)   ;      p_P = p(i)   ;
        mom_L       = rho(i-1)*u(i-1) ;    	rho_L = rho(i-1) ;      u_L = u(i-1) ;      p_L = p(i-1) ;

        lambda_P = abs(u(i)) + sqrt(gamma * p(i) / rho(i));
        lambda_L = abs(u(i - 1)) + sqrt(gamma * p(i - 1) / rho(i - 1));
        lambda_R = abs(u(i + 1)) + sqrt(gamma * p(i + 1) / rho(i + 1));

        lambda_L = max(lambda_L, lambda_P);
        lambda_R = max(lambda_R, lambda_P);
        
        vel_flux_R  = rho_R*u_R*u_R +p_R ;    e_R = e(i+1) ;
        vel_flux_P  = rho_P*u_P*u_P +p_P ;    e_P = e(i)   ;
        vel_flux_L  = rho_L*u_L*u_L +p_L ;    e_L = e(i-1) ;
        
        energy_flux_R   = u_R * ( e_R + p_R ) ;
        energy_flux_P   = u_P * ( e_P + p_P ) ;
        energy_flux_L   = u_L * ( e_L + p_L ) ;
        
        rho_fluxR   = 0.5*( mom_P + mom_R ) -0.5*lambda_R*( rho_R - rho_P ) ;
        rho_fluxL   = 0.5*( mom_P + mom_L ) -0.5*lambda_L*( rho_P - rho_L ) ;
        
        vel_fluxR   = 0.5*( vel_flux_P + vel_flux_R ) -0.5*lambda_R*( mom_R - mom_P ) ;
        vel_fluxL   = 0.5*( vel_flux_P + vel_flux_L ) -0.5*lambda_L*( mom_P - mom_L ) ;
        
        energy_fluxR    = 0.5*( energy_flux_P + energy_flux_R ) -0.5*lambda_R*( e_R - e_P );
        energy_fluxL    = 0.5*( energy_flux_P + energy_flux_L ) -0.5*lambda_L*( e_P - e_L );

        if (i == N - 1)
            if (u(N-1) < a(N-1) && u(N - 1) > 0) 
                vc = u(N - 1) / a(N - 1);

                xd = vc * (1 + gamma) / (2 - vc + vc*gamma);

                uoutlet = u(N - 1) / xd;

                c0 = 0.5 * (1 + gamma) * uoutlet;

                p0 = p(N - 1) * (c0 / a(N - 1)) ^ (2*gamma/ (gamma - 1));

                poutlet = p0 * (2 / (1 + gamma))^(2*gamma/(gamma - 1));

                rhooutlet = gamma * poutlet / uoutlet^2;

                eoutlet = uoutlet.^2*rhooutlet/2 + poutlet / (gamma - 1);

            elseif (u(N - 1) == 0) 
                uoutlet = 2 / (1 + gamma) * a(N - 1);
                poutlet = p(N - 1)*(2 / (1 + gamma)) ^ (2 * gamma / (gamma - 1));
                rhooutlet = gamma * poutlet / uoutlet^2;
                eoutlet = uoutlet.^2*rhooutlet/2 + poutlet / (gamma - 1);
            end

            rho_fluxR = rhooutlet * uoutlet;
            energy_fluxR = uoutlet * (eoutlet + poutlet);
            vel_fluxR = rhooutlet * uoutlet^2 + poutlet;
        
        end
        



        vel_flux    = mom_P - dt/dx(i) * ( vel_fluxR - vel_fluxL ) ;
        
        new_rho(i)  = rho_P - dt/dx(i) * ( rho_fluxR - rho_fluxL ) ;
        
        %vel_flux = mom_P - dt/dx * ( vel_fluxR - vel_fluxL ) - dt*fric + dt*0;
        
        new_u(i)    = vel_flux/new_rho(i) ;
        new_e(i)    = e_P - dt/dx(i) * ( energy_fluxR - energy_fluxL ) ;

        if (i == 2)
            new_rho(1) = rho(1) - dt/dx(i) * rho_fluxL;
            mom1 = rho(1)*u(1) - dt/dx(i) * (vel_fluxL - (0*rho(1)*u(1).^2 +p(1)));
            new_u(1) = mom1 / new_rho(1);
            new_e(1) = e(1) - dt/dx(i) * (energy_fluxL - 0*(u(1) * ( e(1) + p(1) )));
        end
% 
%         if (i == N - 1)
%             new_rho(N) = rho(N) + dt/dx * rho_fluxR;
%             new_u(N) = u(N) + dt/dx * vel_fluxR;
%             new_e(N) = e(N) + dt/dx * vel_fluxR;
%         end

        
        %if (shocktube)
        if (0)
            if (i <= N/2) 
                u_fric = (2 * u_P) / (1+sqrt(1+(16*dt*u_P * lambda) / dia));
        %fric = 0;
            else
                u_fric = u_P;
            end
        else
            %u_fric = (-dia + sqrt(dia*(dia + 16 * dt * u_P* lambda))) / (8 * dt * lambda);
            %fric = rho_P * sign(u_P)*u_P^2 / 8 * lambda / dia * 4.; %dimension Pa/m
            u_fric = (2 * u_P) / (1+sqrt(1+(16*dt*u_P * lambda) / dia));
        end

        %dU_fric = -dt * fric / new_rho(i);
        dU_fric = u_fric - u_P;
        dE_fric = new_rho(i) * ((new_u(i) + dU_fric)^2 - new_u(i)^2) / 2; %kinetic energy lost to friction

        new_u(i) = new_u(i) + dU_fric;
        %new_e(i) = new_e(i) + dE_fric;

        
    end

    end

    if (model_temperature)
        new_T = (new_e - 0.5.*u.*u.*rho) ./ (rho * Cv);
        tdiff = T0 - new_T;

        deltaT = tdiff .* (1 - exp(-dt * wattperkelvinperm3 ./ (rho*Cv)));
        deltaE = deltaT .* rho .* Cv;

        new_e = new_e + deltaE;
    end

    
    new_u(N) = new_u(N - 1);
    
    new_rho(N) = new_rho(N - 1);
    
    new_e(N) = pressureR / (gamma-1) / cfac(end) + 0.5*new_rho(N) .* new_u(N) .* new_u(N);
    
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
