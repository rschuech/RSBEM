function [] = timestepping(input,dump)

timestepping_tic = tic;

switch input.tail.motorBC
    case 'torque'
        fun = @(t,y) yprime_interp_mexed(t,y,dump.best_interpolant,input.performance.nthreads,NaN); %motor_freq isn't known, it's solved for
    case 'freq'
        fun = @(t,y) yprime_interp_mexed(t,y,dump.best_interpolant,input.performance.nthreads,input.tail.motor_freq);
end


clear stop_function
clear global timestepping_convergence
stopfun = @(t,y)stop_function(t,y,dump.y0,input);  %computes convergence criterion on avg speed for when to stop timestepping

[timestepping_solution.x, timestepping_solution.y, timestepping_solution.xe, timestepping_solution.ye, timestepping_solution.ie] = ...
    ode45(fun,[0 input.accuracy.timestepping.Tmax],dump.y0,odeset('InitialStep',input.accuracy.timestepping.initialstep,'events',stopfun,'reltol',input.accuracy.timestepping.reltol,'abstol',input.accuracy.timestepping.abstol));


global timestepping_convergence  % comes from inside ode45, where it was also defined global

plot_timestepping_signal(input,timestepping_convergence);  %update final plot

timings.timestepping = toc(timestepping_tic);
disp([input.paths.namebase.full,'    timestepping took ',num2str(toc(timestepping_tic)/60), ' min','     finished ',datestr(now)]);
save([input.paths.dumpfolder,input.paths.namebase.full,'_timestepping','.mat'],'timestepping_convergence','timestepping_solution','timings');

saveas(input.output.timestepping.fignum,[input.paths.dumpfolder,input.paths.namebase.full,'_timestepping','.fig']);