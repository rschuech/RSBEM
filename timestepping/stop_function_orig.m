function [value,isterminal,direction] = stop_function_orig(t,y,y0,input)

initial_length = 10000;  %initial length of convergence variables
chunk_size = 10000; %how much to grow convergence variables by when necessary (to mitigate slow down due to growing arrays)
maxpts = 5000;
%%

persistent times dists speeds ind

if isempty(dists)  %first time we've entered this function
    dists = NaN(initial_length,1);
    speeds = NaN(initial_length,1);
    times = NaN(initial_length,1);
    ind = 0;
end

ind = ind + 1;   %we have another data point

if ind > length(dists)  %if we're out of space, grow variables
    dists(end+1:end+chunk_size) = NaN;
    speeds(end+1:end+chunk_size) = NaN;
    times(end+1:end+chunk_size) = NaN;
end

dists(ind) = sqrt(sum( (y(1:3) - y0(1:3) ).^2));  %current distance from initial starting point
times(ind) = t;  %current time, duh

if ind > 1 %we have more than one data point
    % p = polyfit(times(1:ind),dists(1:ind),1);  %fitted speed using cumulative data up till now
   
    %non-cumulative appears to converge faster and be less biased than
    %cumulative:
%     p = polyfit(times([1 ind]),dists([1 ind]),1);  %avg speed using only initial and current distances
%     speeds(ind) = p(1);
    speeds(ind) = (dists(ind) - dists(1)) / (times(ind) - times(1));
else
    speeds(ind) = NaN;  %placeholder for first undefined speed value
end

% times_subset = times(max(1, ind - maxpts):ind);
% speeds_subset = speeds(max(1, ind - maxpts):ind);


[times,inds] = sort(times);  % ode45 sometimes goes backward in time, so need to fix data order
%dists = dists(inds);
speeds = speeds(inds);

%following original method doesn't work for speed series that monotonically approach a
%value, which some seem to do....
% diffs = diff(sign(diff(speeds(times <= t))));
% num_mins = sum(diffs == 2); %for maxes, look for -2
% value = min_changes - num_mins;

isterminal = 1; %stop ode45 if avg speed has converged
direction = 0; %stop regardless of convergence from positive or negative direction
%%

pts = speeds(times >= t*(1-input.accuracy.timestepping.datafrac) & times <= t); % take all speed values from last datafrac fraction of data but not after current time (in case ode45 went backwards and there are values after current time)

cv =  nanstd(pts) / nanmean(pts);   %coefficient of variation of above subset of speed points

value = max(input.accuracy.timestepping.minpts - find(times <= t,1,'last'), 0) + (cv - input.accuracy.timestepping.cvtol);  % not converged until both min number of points has been reached and cv is less than cvtol

if rem(ind,input.output.timestepping.plotfreq) == 0
    figure(200)
    plot(times,speeds,'o-');
    hold on
    %plot(times(max(1, ind - maxpts)),speeds(max(1, ind - maxpts)),'ro','markerfacecolor','r','markersize',10);
    hold off
    grid on
    
    title(['cv = ',num2str(cv),'         value = ',num2str(value)]);
    disp(['cv = ',num2str(cv),'         value = ',num2str(value)]);
   % ylim([0.095 0.125]);
    drawnow
    
%     if ~isempty(filename)
%         try
%             print('-dpdf',[filename,'.pdf'])
%         end
%     end

    
end
