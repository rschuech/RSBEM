function [rel_diff] = RMS_fun_diff(fun1, fun2, limits)

% computes the relative difference between two functions (fun2 taken as more accurate than fun1) via the L2 function norm:

% rel_diff = int( sqrt( (fun1 - fun2)^2 ) ) / int( sqrt( fun2 - mean(fun2) )^2 )

% rel_diff = sqrt( int( (fun1 - fun2)^2 ) )  / sqrt( int( fun2 - mean(fun2) )^2 )

% where mean is the average value of a function and both integrals are
% bounded by limits(1) to limits(2)
% 
% fm = @(x) fun2(x);
% fmean = 1/diff(limits) * integral(fm, limits(1),limits(2),'reltol',1E-9,'abstol',1E-12);
% 
% f1 = @(x) sqrt((fun1(x) - fun2(x)).^2);
% num = integral(f1,limits(1),limits(2),'reltol',1E-9,'abstol',1E-12);
% 
% f2 = @(x) sqrt((fun2(x) - fmean).^2);
% den = integral(f2,limits(1),limits(2),'reltol',1E-9,'abstol',1E-12);
% 
% rel_diff = num / den;



fm = @(x) fun2(x);
fmean = 1/diff(limits) * integral(fm, limits(1),limits(2),'reltol',1E-9,'abstol',1E-12);

f1 = @(x) ((fun1(x) - fun2(x)).^2);
num = sqrt(  integral(f1,limits(1),limits(2),'reltol',1E-9,'abstol',1E-12)  );

f2 = @(x) ((fun2(x) - fmean).^2);
den = sqrt(  integral(f2,limits(1),limits(2),'reltol',1E-9,'abstol',1E-12)  );

rel_diff = num / den;