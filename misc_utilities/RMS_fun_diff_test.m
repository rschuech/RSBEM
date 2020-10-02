function [rel_diff] = RMS_fun_diff_test(fun1, fun2, fun_norm, limits)

% computes the relative difference between two functions (fun2 taken as more accurate than fun1) via the L2 function norm:

% used to be
% rel_diff = sqrt( int( (fun1 - fun2)^2 ) )  / sqrt( int( fun_norm - mean(fun_norm) )^2 )
% updated to
% rel_diff = sqrt( int( (fun1 - fun2)^2 ) )  / sqrt( int( fun_norm                  )^2 )
% because we actually do want to normalize by scale of the fun - a
% variation of 1 for something that is around 5 is much more important than
% a variation of 1 for something that is around 1000

% where mean is the average value of a function and both integrals are
% bounded by limits(1) to limits(2)

% orig method uses fun2 as fun_norm, but this weights all kinematics
% components the same when really there are usually some components much
% smaller than others

% so new method normalizes using the (scale of variation of the) vector magnitude of each kinematic
% variable (i.e. U, Omega), making it easier for smaller components to converge (since
% vector magnitude will be much larger than one of the small components)
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



% fm = @(x) fun_norm(x);
% fmean = 1/diff(limits) * integral(fm, limits(1),limits(2),'reltol',1E-9,'abstol',1E-12);

f1 = @(x) ((fun1(x) - fun2(x)).^2);
num = sqrt(  integral(f1,limits(1),limits(2),'reltol',1E-9,'abstol',1E-12)  );


f2 = @(x) ((fun_norm(x) ).^2);
den = sqrt(  integral(f2,limits(1),limits(2),'reltol',1E-9,'abstol',1E-12)  );

rel_diff = num / den;