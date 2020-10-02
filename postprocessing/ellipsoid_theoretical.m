function [fcoeff,frcoeff,tau_a] = ellipsoid_theoretical(a,b,c)

mu = 1E-3 / 1E6; %viscosity in kg / (s micron)

kB = 1.3806488E-23 * (1E6 / 1)^2;  %Boltzmann's constant, kg micron^2 / s^2 / K
T = 298;  %temperature to use for diffusivity calc (K)


Sint = @(x) ( (a.^2 + x).*(b.^2 + x).*(c.^2 + x) ).^(-0.5);

S = integral(Sint,0,Inf);

r = [a; b; c];
Gint = @(x) ( (a.^2 + x).*(b.^2 + x).*(c.^2 + x) ).^(-0.5) .* 1./(r.^2 + x);

G = r.^2 .* integral(Gint,0,Inf,'ArrayValued',true);

vec = 1:3;
for i = vec
    jk = setdiff(vec,i); %other two indices
    H(i,1) = (G(jk(1)) + G(jk(2))) ./ ( r(jk(1)).^2 + r(jk(2)).^2 ) ;% + (S - G(i)) ./  ( r(jk(1)).^2 + r(jk(2)).^2 );
end



fcoeff = 16*pi*mu ./ (S + G); %friction coeff, not force
%F_drag_translation = - fcoeff .* BCs.U  %force = f * speed


frcoeff = 16*pi*mu./3./H; %rotational friction coeff, not torque


tau_a = 1/( kB*T*(1/frcoeff(2) + 1/frcoeff(3)));
%%
%
%
% percenterror_translation = abs(fcoeff_exact - fcoeff) ./ fcoeff_exact * 100;
% percenterror_rotation = abs(frcoeff_exact - frcoeff) ./ frcoeff_exact * 100;
%
% percent_error = [percenterror_translation percenterror_rotation]
