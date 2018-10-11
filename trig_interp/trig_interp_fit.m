function [interpolant] = trig_interp_fit(x,y)
% uses Fourier series to interpolate *evenly spaced* periodic data y(x)
% be sure to leave out final data point that is identical to initial data
% point - this is automatically enforced
% http://www.mathworks.com/help/matlab/math/fast-fourier-transform-fft.html?refresh=true

x = x(:)';  %make sure x is row vector
y = y(:)';  %make sure y is row vector

diffs = diff(x);
deltax = consolidator(diffs,[],[],eps*10);
if ~isscalar(deltax)
    error('x is not evenly spaced')
end

%in case data doesn't start at zero, shift to zero and shift back during
%trig_interp_eval
interpolant.x0 = x(1);
x = x - x(1);

d = fft(y);
m = length(y);  % # data points
interpolant.M = floor((m+1)/2); %order of interpolant

interpolant.a0 = d(1)/m;  %0th order constant
interpolant.an = 2*real(d(2:interpolant.M))/m;  %cos terms
interpolant.bn = -2*imag(d(2:interpolant.M))/m; %sin terms
if rem(m,2) == 0 %even # data points, include final cos term
interpolant.aM = d(interpolant.M+1)/m;
else %odd # data points, leave out final cos term
    interpolant.aM = 0;
end 

interpolant.L = x(end) + deltax;  %period of cyclic data
