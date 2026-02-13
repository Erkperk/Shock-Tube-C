

time_labels = [];

for i=1:6:size(pressure_plot,2)
    %plot(xc(floor(N/2)-10:floor(N/2))-xc(floor(N/2)), pressure_plot(floor(N/2)-10:floor(N/2),i) - 273);
    plot(xc(1:end-1), pressure_plot(1:end-1,i));
%    plot(xc(1:end-1), pressure_plot(1:end-1,i));

    hold on
    %label = sprintf("Time %", floor(pplot_times(i)))
      label= [ "Second " mat2str(floor(pplot_times(i)/1))]

      time_labels = [time_labels; label];
end

legend(time_labels, 'Location', 'northeast')
%legend(time_labels)

xlabel('Distance from opening (m)')
ylabel('Pressure (Pa)')

hold off
