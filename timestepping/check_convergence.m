function [fits] = check_convergence(fits, options)

% if strcmp(mode,'helix')
%     fits.converged.amp = NaN(2,1);
%     fits.converged.R2 = NaN;
% end

switch options.mode
    case 'helix'
        speeds = fits.helix.speed;
    case 'line'
        
        speeds = [fits.line.speed];
end

if ~isfield(options,'direction')
    options.direction = 'first';  %added recently so that we can switch between looking for first time convergence criteria satisfied vs last time
end

diffs = [NaN  abs( diff(speeds) ./ speeds(1:end-1) ) * 100 ];  % percent diffs between consecutive speed fits
good = diffs <= options.diff_cutoff;
good_diff = [false diff(good) == 0];
worked = good_diff & good;
ind = find(worked,1,options.direction);
if ~isempty(ind)
    fits.converged.diff_cutoff = options.diff_cutoff;
    fits.converged.speed = speeds(ind);
    fits.converged.T = fits.T(ind);
    if strcmp(options.mode,'helix')
        fits.converged.R2 = fits.R2(ind);
        if fits.R2(ind) >= min_R2
            fits.converged.amp = fits.helix.amp(:,ind);
        end
    end
    
end



