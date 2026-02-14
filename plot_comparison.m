%% plot_comparison.m — Compare C solver output with Octave workspace
% Run in Octave from the Shock-Tube directory after running the C solver
% and having a saved workspace.mat from the Octave solver.
%
% Usage:  octave --no-gui plot_comparison.m

% ===== Load C solver output =====

% --- Profiles ---
% Format: blocks separated by blank lines, header "# t = <time>"
% Columns: x rho u p T e_int
fid = fopen('profiles.dat', 'r');
c_profiles = {};
c_profile_times = [];
block = [];
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    if length(line) > 2 && line(1) == '#'
        % Header line: extract time
        tok = strsplit(line, '=');
        t_val = str2double(strtrim(tok{2}));
        continue;
    end
    if isempty(strtrim(line))
        % End of block
        if ~isempty(block)
            c_profiles{end+1} = block;
            c_profile_times(end+1) = t_val;
            block = [];
        end
        continue;
    end
    vals = sscanf(line, '%e %e %e %e %f %e');
    if length(vals) == 6
        block = [block; vals'];
    end
end
if ~isempty(block)
    c_profiles{end+1} = block;
    c_profile_times(end+1) = t_val;
end
fclose(fid);

fprintf('Loaded %d profile snapshots from C solver\n', length(c_profiles));

% --- Timeseries ---
% Columns: t massflow u_out rho_out p_out T_out totmass
c_ts = dlmread('timeseries.dat', ' ', 1, 0);  % skip header
fprintf('Loaded %d timeseries rows from C solver\n', size(c_ts, 1));

% ===== Load Octave workspace (if available) =====
have_octave = false;
if exist('workspace.mat', 'file')
    ws = load('workspace.mat');
    have_octave = true;
    fprintf('Loaded Octave workspace\n');
end

% ===== Plot profiles at selected times =====
% Pick ~4 evenly spaced times including near start and end
n_prof = length(c_profile_times);
if n_prof >= 4
    idx_sel = round(linspace(2, n_prof, 4));  % skip t~0
else
    idx_sel = 1:n_prof;
end

figure(1); clf;
for k = 1:length(idx_sel)
    ii = idx_sel(k);
    D = c_profiles{ii};  % columns: x rho u p T e_int

    subplot(2, 2, 1);
    hold on;
    plot(D(:,1)/1000, D(:,4)/1e5, '-', 'DisplayName', sprintf('C t=%.1f', c_profile_times(ii)));
    xlabel('x (km)'); ylabel('Pressure (bar)'); title('Pressure profile');

    subplot(2, 2, 2);
    hold on;
    plot(D(:,1)/1000, D(:,2), '-', 'DisplayName', sprintf('C t=%.1f', c_profile_times(ii)));
    xlabel('x (km)'); ylabel('Density (kg/m^3)'); title('Density profile');

    subplot(2, 2, 3);
    hold on;
    plot(D(:,1)/1000, D(:,3), '-', 'DisplayName', sprintf('C t=%.1f', c_profile_times(ii)));
    xlabel('x (km)'); ylabel('Velocity (m/s)'); title('Velocity profile');

    subplot(2, 2, 4);
    hold on;
    plot(D(:,1)/1000, D(:,5), '-', 'DisplayName', sprintf('C t=%.1f', c_profile_times(ii)));
    xlabel('x (km)'); ylabel('Temperature (K)'); title('Temperature profile');
end

% Overlay Octave results if available
if have_octave && isfield(ws, 'pressure_plot') && isfield(ws, 'xc')
    xc_oct = ws.xc(:);
    n_oct_prof = size(ws.pressure_plot, 2);
    if n_oct_prof >= 4
        oct_sel = round(linspace(2, n_oct_prof, 4));
    else
        oct_sel = 1:n_oct_prof;
    end
    oct_times = ws.pplot_times;

    for k = 1:length(oct_sel)
        jj = oct_sel(k);

        subplot(2, 2, 1); hold on;
        plot(xc_oct/1000, ws.pressure_plot(:, jj)/1e5, '--', 'DisplayName', sprintf('Oct t=%.1f', oct_times(jj)));

        subplot(2, 2, 2); hold on;
        plot(xc_oct/1000, ws.rho_plot(:, jj), '--', 'DisplayName', sprintf('Oct t=%.1f', oct_times(jj)));

        subplot(2, 2, 3); hold on;
        plot(xc_oct/1000, ws.u_plot(:, jj), '--', 'DisplayName', sprintf('Oct t=%.1f', oct_times(jj)));

        subplot(2, 2, 4); hold on;
        plot(xc_oct/1000, ws.T_plot(:, jj), '--', 'DisplayName', sprintf('Oct t=%.1f', oct_times(jj)));
    end
end

for sp = 1:4
    subplot(2, 2, sp);
    legend('Location', 'best');
end

print(gcf, 'comparison_profiles.png', '-dpng', '-r150');
fprintf('Saved comparison_profiles.png\n');

% ===== Plot timeseries =====
figure(2); clf;

subplot(2, 2, 1);
plot(c_ts(:,1), c_ts(:,2), 'b-', 'DisplayName', 'C');
if have_octave && isfield(ws, 'massflow_out')
    hold on;
    plot(ws.plot_times, ws.massflow_out, 'r--', 'DisplayName', 'Octave');
end
xlabel('Time (s)'); ylabel('Mass flow (kg/m^2/s)'); title('Outlet mass flow');
legend('Location', 'best');

subplot(2, 2, 2);
plot(c_ts(:,1), c_ts(:,5)/1e5, 'b-', 'DisplayName', 'C');
if have_octave && isfield(ws, 'p_out')
    hold on;
    plot(ws.plot_times, ws.p_out/1e5, 'r--', 'DisplayName', 'Octave');
end
xlabel('Time (s)'); ylabel('Pressure (bar)'); title('Outlet pressure');
legend('Location', 'best');

subplot(2, 2, 3);
plot(c_ts(:,1), c_ts(:,6), 'b-', 'DisplayName', 'C');
if have_octave && isfield(ws, 'T_out')
    hold on;
    plot(ws.plot_times, ws.T_out, 'r--', 'DisplayName', 'Octave');
end
xlabel('Time (s)'); ylabel('Temperature (K)'); title('Outlet temperature');
legend('Location', 'best');

subplot(2, 2, 4);
plot(c_ts(:,1), c_ts(:,7), 'b-', 'DisplayName', 'C');
if have_octave && isfield(ws, 'totmass')
    hold on;
    plot(ws.plot_times, ws.totmass, 'r--', 'DisplayName', 'Octave');
end
xlabel('Time (s)'); ylabel('Total mass (kg)'); title('Total mass');
legend('Location', 'best');

print(gcf, 'comparison_timeseries.png', '-dpng', '-r150');
fprintf('Saved comparison_timeseries.png\n');
