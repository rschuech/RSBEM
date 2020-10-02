function [fits, options] = path_fitting(t, points, options, fits)

mode = 'line';  %helix mode has been modified and needs to be debugged, tested

%%
if strcmp(mode,'helix')
    min_R2 = 0.98;
    max_iter_helix = 50;
else
    min_R2 = NaN;
end

% bounds_tol_inds_helix =  logical([               1         1        1       1            1           1  1   1             1 1 1]); %untested as of yet


if strcmp(mode,'helix')
    
    %               a1        a2       v       f          phase         shift                    angles
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

objfun = @(x) line_fit(x(1:3), x(4:6), x(7), t, points);
confun = @(x) unit_slope_constraint(x(4:6));

%% line fit

while fits.line.flag(end) <= 0 && fits.line.iters(end) < options.line.max_iter  %fit failed in some way
    fits.line.iters(end) = fits.line.iters(end) + 1;
    
    [options] = choose_guess(options, fits); %updates guess field
    
    try
        [solution_line, rmse_line, fits.line.flag(end)] = fmincon(objfun,options.line.guess,[],[],[],[],options.line.lb,options.line.ub,confun,optimset('tolx',1E-14, 'tolfun', 1E-12, 'tolcon',1E-15,'maxfunevals',5E4, 'maxiter',5E4,'display','off','useparallel',false));
    catch
        solution_line = NaN(7,1);
        rmse_line = NaN;
        fits.line.flag(end) = -20;
    end
    %         fits.line
    
    if any( solution_line(options.line.bound_tol_inds) < options.line.lb(options.line.bound_tol_inds) + abs(options.line.lb(options.line.bound_tol_inds)) * (1 - options.bound_tol))  ||  any( solution_line(options.line.bound_tol_inds) > options.line.ub(options.line.bound_tol_inds) - abs(options.line.ub(options.line.bound_tol_inds)) * (1 - options.bound_tol))
        fits.line.flag(end) = -10;  %means solution is too close to bounds
    end
    
    fits.debug.solution{end}(:,fits.line.iters(end)) = solution_line;
    fits.debug.rmse{end}(fits.line.iters(end)) = rmse_line;
    fits.debug.flag{end}(fits.line.iters(end)) = fits.line.flag(end);
    
    
    %         fits.line.flag = flag;  %record whether fit succeeds or not
    
    if fits.line.flag(end) > 0
        fits.line.intercept(:,end) = solution_line(1:3);
        fits.line.slope(:,end) = solution_line(4:6);
        fits.line.speed(end) = solution_line(7);
    end
    
end



if strcmp(mode,'helix')
    %% helix fit
    
    while fits.R2(c) < min_R2 && fits.helix.iters(c) < max_iter_helix
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
        
        [solution, rmse_helix, flag] = fmincon(objfun,guess,[],[],[],[],helix_lb,helix_ub,[],optimset('tolx',1E-8, 'tolfun', 1E-8, 'tolcon',1E-10,'maxfunevals',3E3, 'maxiter',3E3,'display','off','useparallel',false));
        
        if any( helix_line(bound_tol_inds_helix) < helix_lb(bound_tol_inds_helix) + abs(helix_lb(bound_tol_inds_helix)) * (1 - bound_tol))  ||  any( helix_line(bound_tol_inds_helix) > helix_ub(bound_tol_inds_helix) - abs(helix_ub(bound_tol_inds_helix)) * (1 - bound_tol))
            flag = -10;  %means solution is too close to bounds
        end
        
        %             [rmse_helix] = helix_fit(solution(1), solution(2), solution(3), solution(4), solution(5),solution(6:8), solution(9:11),t, points);
        
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




end

function [options] = choose_guess(options, fits)
if isinf(fits.line.flag(end))  %first try for this T_interrogate
    if length(fits.line.speed) == 1  %very first try ever
        options.line.guess = options.line.guess0;
        %                   guess = line_lb + (line_ub-line_lb).*rand(size(line_lb));
    else %first try for this T, try last solution
        options.line.guess = [fits.line.intercept(:,end-1); fits.line.slope(:,end-1); fits.line.speed(end-1);];
    end
else % line_guess or last solution didn't work, try random guess
    options.line.guess = options.line.lb + (options.line.ub - options.line.lb).*rand(size(options.line.lb));
end

end

