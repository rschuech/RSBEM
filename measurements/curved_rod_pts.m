function [x,y] = curved_rod_pts(u,constants,varargin)


% u = linspace(0,1,500);

r = constants.width / 2;
R = 1 / constants.curvature;
l = constants.length;

den = 2*pi*r + 2*l;

u0 = 0;
u1 = pi*r / den;
u2 = ( pi*r + l/R*(r+R) ) / den;
u3 = ( 2*pi*r + l/R*(r+R) ) / den;
u4 = 1;
% u_cutoffs = [u1 u2 u3];

% u = unique( [u u0 u1 u2 u3 u4]);

t_seg = NaN(size(u));
t_seg( u >= u0 & u <= u1) = 1;
t_seg( u >= u1 & u <= u2) = 2;
t_seg( u >= u2 & u <= u3) = 3;
t_seg( u >= u3 & u <= u4) = 4;

t = NaN(size(u));
t(t_seg == 1) = (2*pi + 2*l/r)*u(t_seg == 1) + pi/2;
t(t_seg == 2) = den/(r+R)*u(t_seg == 2) + 3*pi/2 - pi*r/(r+R);
t(t_seg == 3) = (2*pi + 2*l/r)*u(t_seg == 3) + (-2*l - 3*pi*r)/2/r;
t(t_seg == 4) = den/(r-R)*u(t_seg == 4) + 3*pi/2 - den/(r-R);

x = NaN(size(u));  y = NaN(size(u));
x(t_seg == 1) = r*cos(t(t_seg == 1));                   y(t_seg == 1) = -R + r*sin(t(t_seg == 1));
x(t_seg == 2) = (r+R)*cos(t(t_seg == 2));               y(t_seg == 2) = (r+R)*sin(t(t_seg == 2));
x(t_seg == 3) = R*sin(l/R) + r*cos(t(t_seg == 3));      y(t_seg == 3) = -R*cos(l/R) + r*sin(t(t_seg == 3));
x(t_seg == 4) = (R-r)*cos(t(t_seg == 4));               y(t_seg == 4) = (R-r)*sin(t(t_seg == 4));

% x = x - x(1);
% y = y - y(1);
if length(varargin) == 1
    shift = varargin{1};
    x = x - shift(1);
    y = y - shift(2);
end

% figure(45);  plot(x,y,'-');  hold on
% plot( [ x(u == u0) x(u == u1) ] ,  [ y(u == u0) y(u == u1) ] ,'k-')
% plot( [ x(u == u2) x(u == u3) ] ,  [ y(u == u2) y(u == u3) ] ,'k-')
% hold off
% axis equal
