u = linspace(0,2*pi,500);  h = repmat(2,size(u));
time = 0;
hair_type = 'Coplanar_Hairs';
% hair_type = 'Normal_Top_Hairs';

[~, vel ,~,~] = transverse_hairs_parameterized(u, h, time, geom.transverse , hair_type );

speed = sqrt(sum(vel.^2,1));

figure(493)
histogram(speed,25);