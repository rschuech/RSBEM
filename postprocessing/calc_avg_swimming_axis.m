function [avg_swimming_axis] = calc_avg_swimming_axis(fits, timestepping_solution, Mesh)



[shift,angle,rotvec] = align_path(fits,timestepping_solution.y(1,:));
sol = timestepping_solution.y(timestepping_solution.x <= fits.converged.T , : );
time = timestepping_solution.x(timestepping_solution.x <= fits.converged.T);

[time, inds] = unique( time ) ;  %sometimes the ode45 output contains repeated time points
sol = sol(inds,:);

swimming_axis = NaN(length(time),3);
for i = 1:length(time)
    % this section rotates body orientation to match aligned path
    rotmat = A_1_matrix(sol(i,4:6));  % update to match current time, along orig path
    orientation = rotmat * Mesh(1).orientation;  orientation = orientation ./ repmat(  sqrt(sum(orientation.^2)),3,1);
    rotmat = rotate_arbitrary_vector( rotvec, angle);  % rotate to match aligned path
    orientation = rotmat * orientation;  orientation = orientation ./ repmat(  sqrt(sum(orientation.^2)),3,1);
    
    % now, rotate fixed frame aligned path (x-axis) into current body frame
    rotmat = orientation';  % the coordinate transformation matrix from fixed frame to body frame is simply the transpose of the body frame basis vectors
    
    swimming_axis(i,:) = rotmat * [1 0 0]';   swimming_axis(i,:) =  swimming_axis(i,:)/ sqrt(sum( swimming_axis(i,:).^2));
end

% pp = spline(time,swimming_axis');

 avg_swimming_axis =  ( 1/( time(end) - time(1)) * trapz(time,swimming_axis,1) )'; 

% integral using a spline was taking forever and unclear if it would be any
% more accurate anyway
% avg_swimming_axis = 1/(time(end)-time(1)) * integral( @(t) ppval(pp,t)  , time(1), time(end),'ArrayValued',true);

% renormalize
avg_swimming_axis = avg_swimming_axis / sqrt(sum(avg_swimming_axis.^2));
