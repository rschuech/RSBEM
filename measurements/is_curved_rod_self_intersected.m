function [temp, is_intersected] = is_curved_rod_self_intersected(Length, constants, do_width_correction)

if do_width_correction
    constants.width = 2*Length/(pi - 4)*( sqrt( 1 + (pi - 4)/Length*constants.width ) - 1 );  % from mean width to actual width
end

constants.length = Length - constants.width;
if constants.length < 0
    is_intersected = 1;
    temp = [];
    return
end

[~,~,~,~,~,~,~,shift] = curved_rod_parameters(constants);

%        u = unique( [ linspace(0,1,500)   0 u1 u2 u3 1 ] );
u = linspace(0,1 - 1E-3,500)  ;  %avoid false positive report of intersection at starting/ending point

[x,y] = curved_rod_pts(u,constants,shift);
%      is_intersected = double(~isempty(  InterX([x; y])  )); % self intersections detected
is_intersected = double(~isempty(  intersections(x,y)  )); % self intersections detected

temp = []; % no inequality constraints for fmincon

