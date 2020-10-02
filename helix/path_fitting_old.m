function [fits] = path_fitting_old(timestepping_solution)

mode = 'line';


%%
if strcmp(mode,'helix')
    min_R2 = 0.98;
    max_iter = 50;
else
      min_R2 = NaN;
end



diff_cutoff_1 = 0.2;  %desired cutoff % for convergence
diff_cutoff_2 = 0.5; %if above fails, try this
diff_cutoff_3 = 1; %if above also fails, try this before giving up


clear fits

fits.T = [ 0.03125/8   0.03125/4   0.03125/2   0.03125     0.0625   0.125 0.25  0.5 1 2 4 8 16  32 64 128 ];  %
%shorten T in case timestepping didn't run long enough
while fits.T(end) > timestepping_solution.x(end)
    fits.T(end) = [];
end


%             intercept slope speed
line_guess = [0 0 0   1 0 0  50]';
line_lb = [-50 -50 -50   0 0 0  0]';
line_ub = [50 50 50     1 1 1   500]';


fits.line.intercept = NaN(3,length(fits.T));
fits.line.slope = NaN(3,length(fits.T));
fits.line.speed = NaN(1,length(fits.T));


if strcmp(mode,'helix')
    
    %        a1        a2       v       f     phase           shift                    angles
    helix_lb =    [0          0      0       -2000        -20*1      -10  -10   -10         -pi/4 -pi/4 -pi/4    ]';
    helix_ub =    [10        10     100      2000          20*1        10   10    10          pi/4 pi/4 pi/4      ]';
    
    fits.helix.amp = NaN(2,length(fits.T));
    fits.helix.speed = NaN(1,length(fits.T));
    fits.helix.freq = NaN(1,length(fits.T));
    fits.helix.phase = NaN(1,length(fits.T));
    fits.helix.shift = NaN(3,length(fits.T));
    fits.helix.rotation = NaN(3,length(fits.T));
    
    fits.R2 = -Inf(1,length(fits.T));
    
    fits.helix.iters = zeros(size(fits.R2));
    
end



c = 0;
for T = fits.T
    c = c + 1;
    %     c/length(Ts)
    
    
    if size(timestepping_solution.x,2) == 1
        t = timestepping_solution.x(timestepping_solution.x <= T);
        points = timestepping_solution.y(timestepping_solution.x <= T,1:3);
    else
        t = timestepping_solution.x(timestepping_solution.x <= T)';
        points = timestepping_solution.y(1:3,timestepping_solution.x <= T)';
    end
    
    %% line fit
    
    if c == 1
        guess = line_guess;
    else
        guess = [fits.line.intercept(:,c-1); fits.line.slope(:,c-1); fits.line.speed(c-1);];
    end
    
    objfun = @(x) line_fit(x(1:3), x(4:6), x(7), t, points);
    confun = @(x) unit_slope_constraint(x(4:6));
    
    solution_line = fmincon(objfun,line_guess,[],[],[],[],line_lb,line_ub,confun,optimset('tolx',1E-8, 'tolfun', 1E-8, 'tolcon',1E-10,'maxfunevals',1E6, 'maxiter',1E6,'display','off','useparallel',false));
    
    [rmse_line] = line_fit(solution_line(1:3), solution_line(4:6), solution_line(7), t, points);
    
    fits.line.intercept(:,c) = solution_line(1:3);
    fits.line.slope(:,c) = solution_line(4:6);
    fits.line.speed(c) = solution_line(7);
    
    
    if strcmp(mode,'helix')
        %% helix fit
        
        while fits.R2(c) < min_R2 && fits.helix.iters(c) < max_iter
            fits.helix.iters(c) = fits.helix.iters(c) + 1;
            [T  fits.helix.iters(c)]
            objfun = @(x) helix_fit(x(1),x(2),x(3), x(4), x(5), x(6:8),x(9:11), t,points);
            
            if c == 1 || fits.R2(c-1) < min_R2
                
                guess = helix_lb + (helix_ub-helix_lb).*rand(size(helix_lb));
            else
                %             lo_frac = 0.95;  hi_frac = 1.05;
                %               lo_frac = 0.9;  hi_frac = 1.1;
                lo_frac = 0.75;  hi_frac = 1.25;
                sol = [ fits.helix.amp(:,c-1); fits.helix.speed(:,c-1); fits.helix.freq(:,c-1); fits.helix.phase(:,c-1); fits.helix.shift(:,c-1); fits.helix.rotation(:,c-1); ];
                guess = sol .* (  lo_frac + (hi_frac - lo_frac)*rand(size(helix_lb)) );
            end
            
            solution = fmincon(objfun,guess,[],[],[],[],helix_lb,helix_ub,[],optimset('tolx',1E-8, 'tolfun', 1E-8, 'tolcon',1E-10,'maxfunevals',3E3, 'maxiter',3E3,'display','off','useparallel',false));
            
            [rmse_helix] = helix_fit(solution(1), solution(2), solution(3), solution(4), solution(5),solution(6:8), solution(9:11),t, points);
            
            R2_temp = 1 - rmse_helix^2 / rmse_line^2;
            
            if R2_temp > fits.R2(c)
                fits.R2(c) = R2_temp;
                
                fits.helix.amp(:,c) = solution(1:2);
                fits.helix.speed(:,c) = solution(3);
                fits.helix.freq(:,c) = solution(4);
                fits.helix.phase(:,c) = solution(5);
                fits.helix.shift(:,c) = solution(6:8);
                fits.helix.rotation(:,c) = solution(9:11);
            end
            
        end
        
    end
    
    if c >= 3
        [fits] = check_convergence(fits, diff_cutoff_1,min_R2,mode);
        if ~isnan(fits.converged.speed)
            return
        end
    end
    
    
    
end


% if diff_cutoff_1 failed for all T, try diff_cutoff_2 and _3
[fits] = check_convergence(fits, diff_cutoff_2,min_R2,mode);
if ~isnan(fits.converged.speed)
    return
end

[fits] = check_convergence(fits, diff_cutoff_3,min_R2,mode);
if ~isnan(fits.converged.speed)
    return
end


end




%%
function [fits] = check_convergence(fits,diff_cutoff,min_R2,mode)

fits.converged.diff_cutoff = diff_cutoff;
fits.converged.speed = NaN;
if strcmp(mode,'helix')
    fits.converged.amp = NaN(2,1);
    fits.converged.R2 = NaN;
end

switch mode
    case 'helix'
speeds = fits.helix.speed;
    case 'line'
        speeds = fits.line.speed;
end

diffs = [NaN  abs( diff(speeds) ./ speeds(1:end-1) ) * 100 ];  % percent diffs between consecutive speed fits
good = diffs <= diff_cutoff;
good_diff = [false diff(good) == 0];
worked = good_diff & good;
ind = find(worked,1,'first');
if ~isempty(ind)
    fits.converged.speed = speeds(ind);
    if strcmp(mode,'helix')
        fits.converged.R2 = fits.R2(ind);
        if fits.R2(ind) >= min_R2
            fits.converged.amp = fits.helix.amp(:,ind);
        end
    end
    
end




end

