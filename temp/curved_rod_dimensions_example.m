
AR1 = 1:1:12;
AR2 = 0:0.1:0.9;

[AR1, AR2] = ndgrid(AR1, AR2);  %all combinations of AR1, AR2 above

V = 1;  %micron^3


[minor_radius, major_radius, angle, arclength]  = curved_rod_dimensions(AR1, AR2, V);

% some outputs will be NaN due to impossibility of shapes

NaNs = isnan(minor_radius);  % minor_radius should always be defined so NaN means shape is impossible

minor_radius = minor_radius(~NaNs);
major_radius = major_radius(~NaNs);
angle = angle(~NaNs);
arclength = arclength(~NaNs);

output = [AR1(~NaNs) AR2(~NaNs) minor_radius major_radius angle arclength];

save('curved_rod_dimensions.txt','output','-ascii');  %columns are AR1 AR2 minor_radius major_radius angle arclength