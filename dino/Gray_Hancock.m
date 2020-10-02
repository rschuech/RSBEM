b = 5.9;  %wave amplitude, um
lambda = 23.7;  %wavelength, um
n = 1;  % number of wavelengths
f = 45.7;  % frequency, cycles / s
a = 16.4;  %radius of (spherical) body, um
a = 20
d = 0.15;  %radius of tail, um




V = 2*f*pi^2*b^2/lambda * 1/ ( 1 + 4*pi^2*b^2/lambda^2 - (1+2*pi^2*b^2/lambda^2)^(1/2) * 3*a/n/lambda * ( log(d/2/lambda) + 1) )  %avg. swimming speed, um/s