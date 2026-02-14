%% Plot pressure decline: HLLC vs Lax-Friedrich vs Analytical
% Run both first:  ./shock_tube && ./shock_tube --lf

gamma = 2125/1622;  % Cp/Cv
pL = 165e5;
L  = 4400;

% --- Load numerical data ---
d1 = dlmread('timeseries.dat', ' ', 1, 0);
d2 = dlmread('timeseries_lf.dat', ' ', 1, 0);

% --- Analytical solution (ideal-gas rarefaction reflecting off closed wall) ---
% Estimate real-gas sound speed from wave arrival time
a_L = L / 8.93;   % ~493 m/s
t_arr = L / a_L;

t_an = [0, t_arr - 0.001, t_arr:0.001:20];
p_an = pL * ones(size(t_an));
for k = 1:length(t_an)
    if t_an(k) > t_arr
        a_wall = a_L*(3-gamma)/(gamma+1) + 2*(gamma-1)/(gamma+1) * L/t_an(k);
        a_wall = min(a_wall, a_L);
        p_an(k) = pL * (a_wall / a_L)^(2*gamma/(gamma-1));
    end
end

% --- Plot ---
figure(1); clf;

subplot(2,1,1);
plot(d1(:,1), d1(:,8)/1e5, 'b-', 'LineWidth', 1.5); hold on;
plot(d2(:,1), d2(:,8)/1e5, 'r--', 'LineWidth', 1.5);
plot(t_an, p_an/1e5, 'k-', 'LineWidth', 2);
legend('HLLC (2nd order)', 'Lax-Friedrich (1st order)', 'Analytical (ideal gas)', ...
       'Location', 'southwest');
xlabel('Time (s)'); ylabel('Pressure (bar)');
title('Pressure decline at closed end');
grid on;

subplot(2,1,2);
plot(d1(:,1), d1(:,8)/1e5, 'b-', 'LineWidth', 1.5); hold on;
plot(d2(:,1), d2(:,8)/1e5, 'r--', 'LineWidth', 1.5);
plot(t_an, p_an/1e5, 'k-', 'LineWidth', 2);
legend('HLLC (2nd order)', 'Lax-Friedrich (1st order)', 'Analytical', ...
       'Location', 'southwest');
xlabel('Time (s)'); ylabel('Pressure (bar)');
title('Zoom on arrival');
grid on;
xlim([8.5 10.5]); ylim([130 166]);
