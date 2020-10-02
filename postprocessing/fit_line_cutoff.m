function [fits] = fit_line_cutoff(input,timestepping_solution,T)

%  timestepping_solution = temp.timestepping_solution;
ind = find(timestepping_solution.x > T,1,'first');
timestepping_solution.x = timestepping_solution.x(1:ind);
timestepping_solution.y = timestepping_solution.y(1:ind,:);
timestepping_solution.refpoint = timestepping_solution.refpoint(1:ind,:);

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

  fits.T(end+1) = T;
    fits.line.intercept(:,end+1) = NaN(3,1);
    fits.line.slope(:,end+1) = NaN(3,1);
    fits.line.speed(end+1) = NaN;
    fits.line.iters(end+1) = 0;
    fits.line.flag(end+1) = -Inf;
    
    fits.debug.solution{end+1} = [];
    fits.debug.rmse{end+1} = [];
    fits.debug.flag{end+1} = [];
    
    
    [fits, options] = path_fitting(timestepping_solution.x, timestepping_solution.refpoint, options, fits);
    


