function [] = plot_fits(fits, timestepping_solution)



maxT = 10; %only plot last maxT seconds
minT = 1; %only plot first minT seconds

for i = 1:size(fits.T,2)
    T = fits.T(i);

    if size(timestepping_solution.x,2) == 1
    t = timestepping_solution.x(timestepping_solution.x <= T);
    points = timestepping_solution.y(timestepping_solution.x <= T,1:3);
    else
            t = timestepping_solution.x(timestepping_solution.x <= T)';
        points = timestepping_solution.y(1:3,timestepping_solution.x <= T)';
    end


    inds = t >= t(end) - maxT;
    %       inds = t <= minT;
    %      inds = t >= t(end)*0.5 & t <= t(end)*0.51;
    %      inds = t > 0;


    figure(286)
    p1 =  plot3(points(inds  ,1),points(inds,2),points(inds,3),'.b-');

    [rmse_helix, fitted_helix] = helix_fit(fits.helix.amp(1,i), fits.helix.amp(2,i), fits.helix.speed(:,i), fits.helix.freq(:,i), fits.helix.phase(:,i), fits.helix.shift(:,i), fits.helix.rotation(:,i),t, points);
    [rmse_line, fitted_line] = line_fit(fits.line.intercept(:,i), fits.line.slope(:,i), fits.line.speed(:,i), t, points);

    hold on
    p2 = plot3(fitted_helix(inds,1),fitted_helix(inds,2),fitted_helix(inds,3),'r--','linewidth',1);  axis equal
    p3 = plot3(fitted_line(inds,1),fitted_line(inds,2),fitted_line(inds,3),'b--','linewidth',1);  axis equal

    hold off
    title(['iter = ',num2str(i), '      Tmax = ',num2str(T),'      RMSE = ',num2str(rmse_helix), '     R^2 = ',num2str(fits.R2(i))]);

%     break
    pause

end