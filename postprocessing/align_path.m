function [shift,angle,rotvec] = align_path(fits,y0)
% align swimming path centerline to x-axis


ind = find(fits.line.speed == fits.converged.speed);
intercept = fits.line.intercept(:,ind);
slope = fits.line.slope(:,ind) * fits.line.speed(ind);  % slope and speed are combined in geom3d line representation

% geom3d line object is defined as intercept and slope directly
line = [intercept(:); slope(:)]';


point = projPointOnLine3d(y0(1:3), line);

shift = [0 0 0] - point; %want this point on the centerline to be our new origin

slope_dir = slope / sqrt(sum(slope.^2));
angle = acos( dot(slope_dir,[1 0 0]')  );
rotvec = cross(slope_dir,[1 0 0]');  rotvec = rotvec / sqrt(sum(rotvec.^2));



