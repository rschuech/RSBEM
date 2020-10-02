
function [mean_arclength, arclengths] = tail_arclength(geom, nlambda)

geom.nlambda = nlambda;

% geom.lambda = 22.2;
% geom.nlambda = 1.87;
% geom.amp = [1.752  11.4/2];
% geom.kE = [3 * 0.24    2 * 0.05 * 1.5];
% geom.omega = 2 * pi / 0.0219 ;
% 
% % geom.t_transition =  3.20735258036;  %thick tail angled
% geom.t_transition = 3.05064508356;  %from revised geom
% 
% geom.translation = ([ 6.36348374047 3.5527136788e-15 15.2])';  %thin tail angled
% geom.translation = ([ 6.04563244435 3.5527136788e-15 15.2])';   %thick tail angled
% % geom.translation = ([ 1.04563244435 3.5527136788e-15 87.2])';   %thick tail angled
% 
% geom.rotation = [0 0 pi]';
% geom.radius = 0.15   * 1;
 geom.t_max = 2*pi * geom.nlambda;  %t value where hemispherical end cap is centered



% geom.tail_angle = -25 *0 ;
% 
% geom.rotation_pt = [6.85279653889 3.5527136788e-15 15.2]';  %thin tail angled
% geom.rotation_pt = [6.85279653889 3.5527136788e-15 15.2]';  %thick tail angled
% 
% geom.rotation_vec = [0 0 1];  %always vertically up along z dir


%%
c = 0;  clear arclengths
T = linspace(0,10,50);
parfor c = 1:length(T)
    time = T(c);
    integrand_fun = @(t) integrand_wrapper(t, time, geom);
    
    %      fun = @(t2) arrayfun(@(t)
    
    arclengths(c) = integral(integrand_fun,0,geom.t_max,'abstol',1E-10,'reltol',1E-6);
end

% f = @(x,y)exp(-hypot(x,y));
% lim1 = -inf;
% lim2 = inf;
% g = @(y1)arrayfun(@(y)integral(@(x)f(x,y),lim1,lim2),y1);
% y = linspace(-10,10);
% plot(y,g(y));
mean_arclength = mean(arclengths);