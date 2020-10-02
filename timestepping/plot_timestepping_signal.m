function [] = plot_timestepping_signal(input,timestepping_convergence)

plot_every_nth = 10;  %only plot every nth datapoint

figure(input.output.timestepping.fignum)
% plot(timestepping_convergence.times(ceil(cutoff_ind*plotfrac):end),timestepping_convergence.speeds(ceil(cutoff_ind*plotfrac):end),'o-');
plot(timestepping_convergence.times(1:plot_every_nth:end),timestepping_convergence.speeds(1:plot_every_nth:end),'o-');
if timestepping_convergence.t > timestepping_convergence.interrogation_interval
    hold on
    plot(timestepping_convergence.times(timestepping_convergence.cutoff_ind),timestepping_convergence.speeds(timestepping_convergence.cutoff_ind),'ro','markerfacecolor','r','markersize',10);
    hold off
end
grid on

title(['error estimate = ',num2str(timestepping_convergence.error_estimate),'   ','freq = ',num2str(timestepping_convergence.motor_freq), '   avg. speed = ',num2str(timestepping_convergence.avg_speed)]);

disp(['t = ',num2str(timestepping_convergence.t),'         error estimate = ',num2str(timestepping_convergence.error_estimate)]);


drawnow

if input.output.timestepping.doprint
    try
        print('-dpdf',[input.paths.dumpfolder,input.paths.namebase.full,'_timestepping','.pdf'])
    catch
        disp(['Couldn''t write to ',[input.paths.dumpfolder,input.paths.namebase.full,'_timestepping','.pdf']]);
    end
end
