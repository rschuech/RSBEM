function [r,R,l,den,u1,u2,u3,shift, midpt] = curved_rod_parameters(constants)

r = constants.width / 2;
R = 1 / constants.curvature;
l = constants.length;

den = 2*pi*r + 2*l;


u1 = pi*r / den;
u2 = ( pi*r + l/R*(r+R) ) / den;
u3 = ( 2*pi*r + l/R*(r+R) ) / den;

[x,y] = curved_rod_pts( (0+u1)/2,constants);
shift = [x y]; % attempt to fix finite precision problem (straight rods located at huge y coords) by shifting geometry close to origin


midpt = [ R*cos(-pi/2 + l/R/2)    ,     R*sin(-pi/2 + l/R/2) ]  - shift  ;  % mid point of centerline