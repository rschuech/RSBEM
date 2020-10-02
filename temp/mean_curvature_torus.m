a = 2;  
c = 4;

v = linspace(0,2*pi,1E6);

mean_curv = (  - cos(v)./(c + a .* cos(v))  + 1/a   ) / 2;
abs_mean_curv = abs(   cos(v)./(c + a .* cos(v))  + 1/a   ) / 2;

% mean(mean_curv)

figure(394);

plot(v,abs_mean_curv,'o-');
grid on
% [  max(mean_curv)   min(mean_curv)  ]

% max(abs_mean_curv)

% abs([ (2*a - c)/2/a/(a-c)      (2*a+c)/2/a/(a+c)  ])  % these are putative extreme values gotten by making cos(v) = -1 or +1 in mean_curv expression above


k1 = cos(v)./(c+a*cos(v));
mean(k1)

figure(3843)

plot(v,k1,'o-')
grid on