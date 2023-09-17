plot(plot_times(2:end)/60,-diff(totmass)./diff(plot_times))
xlabel('Time (min)')
ylabel('Mass flow (kg/s)')
grid minor
