geom0 = [ [Results0.AR1]' [Results0.AR2]'  [Results0.amp]'  [Results0.lambda]'  [Results0.nlambda]'  ];

geom = [ [Results.AR1]' [Results.AR2]'  [Results.amp]'  [Results.lambda]'  [Results.nlambda]'  ];


[shat ,ia,ib]= intersect(geom,geom0,'rows')


%%

nth = 1E5;  nth =8;  nth = 2;
for ind = length(fits.T) - 1
    nth = 100  ;
t1 = 0000;  t2 = fits.T(ind); %t2 = 65536;
t = timestepping_solution.x(timestepping_solution.x >= t1 & timestepping_solution.x <= t2);
x = timestepping_solution.y(timestepping_solution.x >= t1 & timestepping_solution.x <= t2,1);
y = timestepping_solution.y(timestepping_solution.x >= t1 & timestepping_solution.x <= t2,2);
z = timestepping_solution.y(timestepping_solution.x >= t1 & timestepping_solution.x <= t2,3);
t = t(1:nth:end);  x = x(1:nth:end); y = y(1:nth:end);  z = z(1:nth:end);

[~,fitted] = line_fit(fits.line.intercept(:,ind),fits.line.slope(:,ind),fits.line.speed(ind),t,NaN(length(t),3));   
figure(ind+820)

clf
hold off
plot3(x,y,z,'o-');  
hold on; 
plot3(fitted(:,1),fitted(:,2),fitted(:,3),'-or');   
grid on
title(['ind = ',num2str(ind)]);
drawnow
pause
end