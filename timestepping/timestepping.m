function [fits, timestepping_solution, timings] = timestepping(input,dump,timestepping_solution0,refpoint0,Mesh)

if isfield(input.paths,'namebase')
    name = input.paths.namebase.full;
else
    name = input.paths.fullnamebase;
end

if ~isempty(timestepping_solution0)
    if size(timestepping_solution0.y,1) < size(timestepping_solution0.y,2) %row is variable, make row timestep
        timestepping_solution0.y = timestepping_solution0.y';
    end
end

timestepping_tic = tic;

switch input.bugtype
    case 'bacteria'
        switch input.tail.motorBC
            case 'torque'
                fun = @(t,y) yprime_interp(t,y,dump.best_interpolant,input.kinematics_interp_method,input.periodic_kinematics,NaN); %motor_freq isn't known, it's solved for
            case 'freq'
                fun = @(t,y) yprime_interp(t,y,dump.best_interpolant,input.kinematics_interp_method,input.periodic_kinematics,input.tail.motor_freq);
        end
    case 'dino'
%         temp.x0 = NaN; temp.M = NaN; temp.a0 = NaN; temp.an = NaN; temp.bn = NaN; temp.aM = NaN; temp.L = NaN;
%         fun = @(t,y)
%         yprime_interp_mex(t,y,temp,dump.best_interpolant,input.kinematics_interp_method,input.phase_speed);
%         % can't mex anymore due to use of spline() instead of Fourier
%         interp?
 fun = @(t,y) yprime_interp(t,y,dump.best_interpolant,input.kinematics_interp_method,input.periodic_kinematics,input.phase_speed);
        
end
options.bound_tol = 0.95;  %if solution contains one or more components beyond this fraction of a lower or upper bound, flag the solution as probably bad - need to adjust bounds for the future
options.line.bound_tol_inds = logical([1 1 1 0 0 0 1]);  %only check these solution parameters for being close to bounds (slope of line is exempt from check since likely to be near 1 for x-axis)

switch input.bugtype
    case 'bacteria'
%         options.line.guess0 = [0      0   0    1  0  0   12]';
%         options.line.lb =    [-5000 -10 -10   -1 -1 -1    0]';
%         options.line.ub =    [5000   10  10    1  1  1  100]';
        options.line.guess0 = [Mesh(2).refpoints(:,1)'    Mesh(2).orientation(:,1)'     12]';
        options.line.lb =    [-5000 -5000 -5000   -1 -1 -1    0]';
        options.line.ub =    [5000   5000  5000    1  1  1  100]';
    case 'dino'
%           options.line.guess0 = [0      0   0    1  0  0   150]';
%         options.line.lb =    [-5000 -1000 -1000   -1 -1 -1    0]';
%         options.line.ub =    [5000   1000  1000    1  1  1  1000]';
        options.line.guess0 = [0      0   0    1  0  0   150]';
        options.line.lb =    [-5000 -5000 -5000   -1 -1 -1    0]';
        options.line.ub =    [5000   5000  5000    1  1  1  1000]';
end

options.diff_cutoff = input.accuracy.timestepping.diff_tols(1);  %desired cutoff % for convergence
% diff_cutoff_2 = 0.5; %if above fails, try this
% diff_cutoff_3 = 1; %if above also fails, try this before giving up
options.mode = 'line';
options.line.max_iter = 10;


fits.T = [];
fits.line.intercept = [];
fits.line.slope = [];
fits.line.speed = [];
fits.line.iters = [];
fits.line.flag = [];

fits.converged.diff_cutoff = NaN;
fits.converged.speed = NaN;

fits.debug.solution = {};
fits.debug.rmse = {};
fits.debug.flag = {};


timestepping_solution.x = []; %time
timestepping_solution.y = []; % translation, rotation, and possibly tail phase
timestepping_solution.refpoint = []; %tracked position of refpoint, which is different than translation
converged = false;

for ti = 1:length(input.accuracy.timestepping.T_interrogate)
    if input.performance.verbose
        disp(['Began T_interrogate = ',num2str(input.accuracy.timestepping.T_interrogate(ti)),' at ',datestr(now)]);
    end
    %initialize shat for current T_interrogate
    fits.T(end+1) = input.accuracy.timestepping.T_interrogate(ti);
    fits.line.intercept(:,end+1) = NaN(3,1);
    fits.line.slope(:,end+1) = NaN(3,1);
    fits.line.speed(end+1) = NaN;
    fits.line.iters(end+1) = 0;
    fits.line.flag(end+1) = -Inf;
    
    fits.debug.solution{end+1} = [];
    fits.debug.rmse{end+1} = [];
    fits.debug.flag{end+1} = [];
    
    if ~isempty(timestepping_solution0) && timestepping_solution0.x(end) >= input.accuracy.timestepping.T_interrogate(ti)
        % this doesn't seem to be used anymore?
        error('old code reached')
        timestepping_solution.x = timestepping_solution0.x(timestepping_solution0.x <= input.accuracy.timestepping.T_interrogate(ti));
        timestepping_solution.y = timestepping_solution0.y(timestepping_solution0.x <= input.accuracy.timestepping.T_interrogate(ti),:);
       
    else
        
        if ti == 1
            
            [x, y] = ode45(fun,[0 input.accuracy.timestepping.T_interrogate(ti)],dump.y0,odeset('InitialStep',input.accuracy.timestepping.initialstep,'reltol',input.accuracy.timestepping.reltol,'abstol',input.accuracy.timestepping.abstol));
            
        else
            
            if ~isempty(timestepping_solution0) && timestepping_solution0.x(end) >  timestepping_solution.x(end)
                % this doesn't seem to be used anymore?
                error('old code reached')
                timestepping_solution.x = timestepping_solution0.x;
                timestepping_solution.y = timestepping_solution0.y;
            end
            
            npoints = 100;  safety_factor = 10;
            if length(timestepping_solution.x) > npoints
                step0 = mean(  diff(  timestepping_solution.x(end-100 : end) ) ) / safety_factor;
            else
                step0 = mean(  diff(  timestepping_solution.x ) ) / safety_factor;
            end
         
            %   disp(['computing t = ',num2str(input.accuracy.timestepping.T_interrogate(ti-1)),'     to      t = ',num2str(input.accuracy.timestepping.T_interrogate(ti))]);
            [x, y] = ode45(fun,[timestepping_solution.x(end) input.accuracy.timestepping.T_interrogate(ti)], timestepping_solution.y(end,:)' ,odeset('InitialStep',step0,'reltol',input.accuracy.timestepping.reltol,'abstol',input.accuracy.timestepping.abstol));
            
        end
        
        timestepping_solution.x = [timestepping_solution.x; x];
        timestepping_solution.y = [timestepping_solution.y; y];
        
        refpoints = NaN(length(x),3);
        for ri = 1:length(x)
            refpoints(ri,:) = A_1_matrix(y(ri,4:6)) * refpoint0   +  y(ri,1:3)';
        end
        timestepping_solution.refpoint = [timestepping_solution.refpoint; refpoints];
        
    end
    
    [fits, options] = path_fitting(timestepping_solution.x, timestepping_solution.refpoint, options, fits);
    
    if ti >= 3  %need to compare at least two consecutive differences so need at least 3 data points
        [fits] = check_convergence(fits, options);
        %         fits.line.speed
        if ~isnan(fits.converged.speed)
            converged = true;
            break
        else
            save([input.paths.dumpfolder,name,'_partial','.mat'], 'timestepping_solution', 'fits');
        end
    end
    
end

if ~converged % try ever larger convergence tols until it works
    for i = 1:length(input.accuracy.timestepping.diff_tols)
        options.diff_cutoff = input.accuracy.timestepping.diff_tols(i);  options.direction = 'last';
        [fits] = check_convergence(fits, options);
        if ~isnan(fits.converged.speed)
            converged = true;
            %             t_goodenough(f) = temp2.fits.T( temp2.fits.line.speed == temp2.fits.converged.speed);
            %             cutoff_goodenough(f) = speed_tols(i);
            break
        end
    end
end


timings.timestepping = toc(timestepping_tic);
if input.performance.verbose
    disp([name,'    timestepping took ',num2str(toc(timestepping_tic)/60), ' min','     finished ',datestr(now)]);
end
save([input.paths.dumpfolder,name,'_timestepping','.mat'], 'timestepping_solution', 'fits', 'timings');

if exist([input.paths.dumpfolder,name,'_partial','.mat'])
    delete([input.paths.dumpfolder,name,'_partial','.mat']);
end

% saveas(input.output.timestepping.fignum,[input.paths.dumpfolder,input.paths.namebase.full,'_timestepping','.fig']);


%%
