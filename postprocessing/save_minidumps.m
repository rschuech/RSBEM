folder = 'C:\Users\rudi\Desktop\RD\swept_dumps\';
folder2 = 'C:\Users\rudi\Desktop\RD\minidumps\';

files = dir([folder,'*timestepping.mat']);
files = {files.name};

for f = 1:length(files)
    file = files{f};
    f/length(files)
    
    newfile = [file(1:end-16) 'converged_speed.mat'];
%     
    if exist([folder2,newfile],'file')
        continue
    end
    try
    temp = load([folder,file],'timestepping_convergence');
    catch
        continue
    end
    
    
    
    timestepping_convergence = temp.timestepping_convergence;
    timestepping_convergence = rmfield(timestepping_convergence,{'speeds','times'});
    
    %newfile = [file(1:end-16) 'converged_speed.mat'];
    %% compute time-averaged omega over a rotation cycle (and thus for entire trajectory, since omega is in the body frame and thus periodic)
    interp_file = [ file(1:end-16) 'dump.mat'];  %remove "timestepping" and replace with "dump" for interpolation file
    temp = load([folder,interp_file],'interpolant');
    interpolant = temp.interpolant(end).omega;  %last and most accurate interpolant for tail rotation rate omega
    
    fun = @(t,theta) trig_interp_eval(interpolant,theta);
    
    stopfun = @(t,theta) stop_function_theta(t,theta);
    
    sol = ode45(fun,[0 10],0,odeset('InitialStep',1E-12,'events',stopfun,'reltol',1E-10,'abstol',1E-12));
    
    figure(383)
    plot(sol.x,sol.y,'o-','markerfacecolor','b'); grid on;  drawnow
    
    T = sol.xe;  %temporal period of a tail rotation (we don't know this a priori since the BC was constant torque, with variable omega)
    timestepping_convergence.avg_omega = 2*pi / T; %avg omega is simply 2*pi radians / time it took to rotate that much
    %%
    
    
    
    save([folder2,newfile],'timestepping_convergence');
    
end