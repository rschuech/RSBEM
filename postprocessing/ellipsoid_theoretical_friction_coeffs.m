function [f_translation, f_rotation] = ellipsoid_theoretical_friction_coeffs(geom, mu)


% geom.a = 2;    % ellipsoid radius 1 (um)
% geom.b = 0.5;  % ellipsoid radius 2 (um)
% geom.c = 0.5;  % ellipsoid radius 3 (um)  (for a rotationally symmetric ellipsoid, geom.b = geom.c)
% mu = 1E-9;  % dynamic viscosity of water  (  kg / (um s)  )



% following is from Dusenbery appendix

Sint = @(x) ( (geom.a.^2 + x).*(geom.b.^2 + x).*(geom.c.^2 + x) ).^(-0.5);  % (1/um^3)

S = integral(Sint,0,Inf);  % units of integral are units of Sint * x [which must be in um^2] so S is (1/um)

r = [geom.a; geom.b; geom.c];
Gint = @(x) ( (geom.a.^2 + x).*(geom.b.^2 + x).*(geom.c.^2 + x) ).^(-0.5) .* 1./(r.^2 + x);

G = r.^2 .* integral(Gint,0,Inf,'ArrayValued',true);

vec = 1:3;
for i = vec
    jk = setdiff(vec,i); %other two indices
    H(i,1) = (G(jk(1)) + G(jk(2))) ./ ( r(jk(1)).^2 + r(jk(2)).^2 ) ;% + (S - G(i)) ./  ( r(jk(1)).^2 + r(jk(2)).^2 );
end



f_translation = 16*pi*mu ./ (S + G); %friction coeffs (3 x 1 vector) (kg/s) for translation along direction a, b, and c

% f_rand = 3 / sum( 1./fcoeff_exact )  %single effective friction coeff (kg/s) for avg randomized orientation due to rotational Brownian motion

% f_rand (kg/s) * speed (um/s) = drag force (uN)



 f_rotation = 16*pi*mu./3./H; %rotational friction coeffs

