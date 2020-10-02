function [arclength] = bacterial_tail_arclength(tail)

arclength = NaN(size(tail,1),1);
for i = 1:size(tail,1)
    if any(isnan(tail(i,:)))
        continue
    end
    
    
amp = tail(i,1);  lambda = tail(i,2);  nlambda = tail(i,3);
% http://mathworld.wolfram.com/ArcLength.html
% derivatives of parametric flagellum equation from Maple

k = 2*pi./lambda;  KE = k; % as per Shum et al
integrand = @(xi) sqrt( (1)^2     +    (2*amp.*KE.^2.*xi.*exp(-KE.^2.*xi.^2).*cos(k.*xi)-amp.*(1-exp(-KE.^2.*xi.^2)).*sin(k.*xi).*k).^2    +    (2*amp.*KE.^2.*xi.*exp(-KE.^2.*xi.^2).*sin(k.*xi)+amp.*(1-exp(-KE.^2.*xi.^2)).*cos(k.*xi).*k).^2   );
xi_max = nlambda * lambda;
arclength(i) = integral(integrand,0,xi_max);

end