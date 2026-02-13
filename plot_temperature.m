

time_labels = [];

for i=1:20:size(T_plot,2)
    %plot(xc(floor(N/2)-10:floor(N/2))-xc(floor(N/2)), T_plot(floor(N/2)-10:floor(N/2),i) - 273);
    plot(xc(end-100:end-1), T_plot(end-100:end-1,i) - 273);

    hold on
    %label = sprintf("Time %", floor(pplot_times(i)))
      label= [ "Minute " mat2str(floor(pplot_times(i)/60))]

      time_labels = [time_labels; label];
end

legend(time_labels, 'Location', 'southwest')
%legend(time_labels)

xlabel('Distance from opening (m)')
ylabel('Temperature (C)')

hold off
