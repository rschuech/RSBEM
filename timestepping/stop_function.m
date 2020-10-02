function [value,isterminal,direction] = stop_function(t,y,y0,input)

skip_stop_check = true;

safety_factor = 1.1;  %leave this much more data than exactly interrogation_interval when cutting off old data - prevents NaN error_estimate being erroneously reported

% initial_length = 10000;  %initial length of convergence variables
% chunk_size = 10000; %how much to grow convergence variables by when necessary (to mitigate slow down due to growing arrays)

isterminal = 1; %stop ode45 if avg speed has converged
direction = 0; %stop regardless of convergence from positive or negative direction
%%

global timestepping_convergence
persistent ind last_phase_angle plot_ind

if isempty(timestepping_convergence)  %first time we've entered this function, initialize global/persistent variables
    timestepping_convergence.speeds = NaN(input.performance.timestepping.initial_length,1);
    timestepping_convergence.times = NaN(input.performance.timestepping.initial_length,1);
    timestepping_convergence.avg_speed = NaN;
    ind = 0;
    plot_ind = 0;
end


switch input.tail.motorBC
    case 'torque'  %need to compute an estimate of motor freq
        if ind > 0
            timestepping_convergence.motor_freq = ( y(7) - last_phase_angle ) / ( t - timestepping_convergence.times(ind) );
            timestepping_convergence.interrogation_interval = input.accuracy.timestepping.normalized_interrogation_interval ./ timestepping_convergence.motor_freq;
        else
            timestepping_convergence.interrogation_interval = Inf; %ensures we return without hitting the stop condition
        end
        last_phase_angle = y(7);
    case 'freq' %we already know motor freq, it's constant
        timestepping_convergence.motor_freq = input.tail.motor_freq;
        timestepping_convergence.interrogation_interval = input.accuracy.timestepping.normalized_interrogation_interval / timestepping_convergence.motor_freq;
        
end



ind = ind + 1;   %we have another data point

if ind > length(timestepping_convergence.speeds)  %because stored data can get pretty huge, first try deleting old data, then grow vectors once in a while if still needed
    timestepping_convergence.cut_ind = find(timestepping_convergence.times < (t - timestepping_convergence.interrogation_interval*safety_factor), 1, 'last');  %everything before this is pretty old so try shrinking vectors first
    if ~isempty(timestepping_convergence.cut_ind)
        timestepping_convergence.times(1:timestepping_convergence.cut_ind) = [];    timestepping_convergence.speeds(1:timestepping_convergence.cut_ind) = [];  %cut off old data, shrinking vectors
        ind = ind - timestepping_convergence.cut_ind;  %need to shift current ind down after cutting off old data
    end
    
    if ind > length(timestepping_convergence.speeds) %if after cutting off old data we still don't have enough space, grow vectors
        timestepping_convergence.speeds(end+1:end + input.performance.timestepping.chunk_size) = NaN;
        timestepping_convergence.times(end+1:end + input.performance.timestepping.chunk_size) = NaN;
    end
    
end

timestepping_convergence.times(ind) = t;  %current time, duh
timestepping_convergence.t = t;

if ind > 1 && ~skip_stop_check %we have more than one data point and we're not skipping stop checks
    % p = polyfit(times(1:ind),dists(1:ind),1);  %fitted speed using cumulative data up till now
    
    %non-cumulative approach appears to converge faster and be less biased than
    %cumulative:
    %     p = polyfit(times([1 ind]),dists([1 ind]),1);  %avg speed using only initial and current distances
    %     speeds(ind) = p(1);
    dist = sqrt(sum( (y(1:3) - y0(1:3)).^2)); %current distance from initial starting point
    timestepping_convergence.speeds(ind) = (dist ) / (t);  %assumes dist(t = 0) = 0, which it certainly should, and also that simulation starts at t = 0
else
    timestepping_convergence.speeds(ind) = NaN;  %placeholder for first undefined speed value
end


% so far, all timeseries have eventually become decaying periodic signals.
% if this doesn't happen, first try refining ode45 reltol, since this could
% be inaccurate enough to produce a totally different signal
% also, a somewhat too big ode 45 reltol will often cause translational shifts or phase errors in the
% timeseries - can compare unrefined to refined and look for phase
% differences


%%
if t <= timestepping_convergence.interrogation_interval || skip_stop_check %don't have enough data yet, keep going
    value = realmax;  %using Inf sometimes causes events checker to fuck up cause it tries to do 0 * this which is NaN for 0 * Inf
    timestepping_convergence.error_estimate = realmax;
    timestepping_convergence.cutoff_ind = 1; %for plotting purposes
else
    timestepping_convergence.cutoff_time = t - timestepping_convergence.interrogation_interval;
    timestepping_convergence.cutoff_ind = find(timestepping_convergence.times < timestepping_convergence.cutoff_time, 1, 'last');
    
    speeds_subset = timestepping_convergence.speeds(timestepping_convergence.cutoff_ind:ind);
    
    range = [min(speeds_subset)  max(speeds_subset)];
    if ~isempty(range)
        timestepping_convergence.error_estimate = (range(2) - range(1)) / mean(speeds_subset) / 2;  %like coeff of variation except uses the max range of the data instead of stdev
        % / 2 because the estimate is really way too conservative otherwise,
        % since true value should be very close to (max + min)/2 due to
        % decaying periodic nature of the timeseries
        
        timestepping_convergence.avg_speed = (range(1) + range(2)) / 2;  %current best estimate for true avg speed
    else
        timestepping_convergence.error_estimate = realmax;
        timestepping_convergence.avg_speed = NaN;
    end
    
    value = timestepping_convergence.error_estimate - input.accuracy.timestepping.error_tol;
    
    
    
end


plot_ind = plot_ind + 1;

if rem(plot_ind,input.output.timestepping.plotfreq) == 0  %only plot once in a while
    
    
%     timestepping_convergence
%     error_estimate
%     t
%     cutoff_ind
    
    plot_timestepping_signal(input,timestepping_convergence);
    
    
end
