function y = trig_interp_eval(interpolant, x)
% evaluates trigonometric interpolant at input x values (no restriction on
% spacing)
% note that aM = 0 for odd # original data points
% http://www.mathworks.com/help/matlab/math/fast-fourier-transform-fft.html?refresh=true

x = x(:)';  %make sure x is row vector

shifted_x = x - interpolant.x0;  %interpolant is defined for data starting at zero, so shift new x to match

n = 1:length(interpolant.an);
y = interpolant.a0 + interpolant.an*cos(2*pi*n'*shifted_x/interpolant.L) ...
    + interpolant.bn*sin(2*pi*n'*shifted_x/interpolant.L) + interpolant.aM*cos(2*pi*interpolant.M*shifted_x/interpolant.L);
