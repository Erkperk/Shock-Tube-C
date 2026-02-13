if ~exist('plot_times', 'var'), load('workspace.mat'); end
plot((plot_times)/(1), p_base / 100000)   %+61440 NS1, 7380NS1 LT
xlabel('Time (seconds)')
ylabel('Pressure (bar)')
grid on
grid minor
