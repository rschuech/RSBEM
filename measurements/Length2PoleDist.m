function [PoleDist] = Length2PoleDist(Length, constants, check_intersections, do_width_correction)

if do_width_correction
    constants.width = 2*Length/(pi - 4)*( sqrt( 1 + (pi - 4)/Length*constants.width ) - 1 );  % from mean width to actual width
end

constants.length = Length - constants.width;


[~,~,~,~,u1,u2,u3,shift] = curved_rod_parameters(constants);

if check_intersections
    u = linspace(0,1 - 1E-3,500);
    [x,y] = curved_rod_pts(u,constants,shift);
    if ~isempty(  intersections(x,y)  ) % self intersections detected
        PoleDist = NaN;
        return
    end
end
%%

u_poles = [u1/2  (u2+u3)/2];  % the poles

[x,y] = curved_rod_pts(u_poles,constants,shift);
PoleDist = sqrt( diff(x,[],2).^2 + diff(y,[],2).^2 );
