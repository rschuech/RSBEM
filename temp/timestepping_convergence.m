
filter = timestepping_solution.x < 8 & timestepping_solution.x > 0 ;

fun = @(s) mod(s+2*pi,2*pi);
% fun = @(s) s*180/pi;
% fun = @(s) s;

delta = 1;

temp = false(size(filter));  temp(1:delta:end) = true;
filter = filter & temp;

figure(458)
clf
subplot(4,1,1)
plot(timestepping_solution.x(filter),fun(timestepping_solution.y(filter,4)),'o-');
% ylim([0 2*pi]);
ax1 = gca;
grid on
subplot(4,1,2)
plot(timestepping_solution.x(filter),fun(timestepping_solution.y(filter,5)),'o-');
% ylim([0 2*pi]);
ax2 = gca;
grid on
subplot(4,1,3)
plot(timestepping_solution.x(filter),fun(timestepping_solution.y(filter,6)),'o-');
% ylim([0 2*pi]);
ax3 = gca;
grid on

subplot(4,1,4)
plot(timestepping_solution.x(filter),fun(timestepping_solution.y(filter,7)),'o-');
% ylim([0 2*pi]);
ax4 = gca;
grid on
linkaxes([ax1 ax2 ax3 ax4],'x');

%%
tol = 2E-2;
% ind = find( timestepping_solution.x > .1 & ...
%     ( mod(timestepping_solution.y(:,4)+2*pi,2*pi) < tol | mod(timestepping_solution.y(:,4)+2*pi,2*pi) > (2*pi - tol) ) & ...
%     ( mod(timestepping_solution.y(:,5)+2*pi,2*pi) < tol | mod(timestepping_solution.y(:,5)+2*pi,2*pi) > (2*pi - tol) ) & ...
%     ( mod(timestepping_solution.y(:,6)+2*pi,2*pi) < tol | mod(timestepping_solution.y(:,6)+2*pi,2*pi) > (2*pi - tol) )        ,1,'first');

ind = find( timestepping_solution.x > .1 & ...
    ( mod(timestepping_solution.y(:,4)+2*pi,2*pi) < tol | mod(timestepping_solution.y(:,4)+2*pi,2*pi) > (2*pi - tol) ) & ...
    ( mod(timestepping_solution.y(:,5)+2*pi,2*pi) < tol | mod(timestepping_solution.y(:,5)+2*pi,2*pi) > (2*pi - tol) ) & ...
    ( mod(timestepping_solution.y(:,6)+2*pi,2*pi) < tol | mod(timestepping_solution.y(:,6)+2*pi,2*pi) > (2*pi - tol) ) & ...
    ( mod(timestepping_solution.y(:,7)+2*pi,2*pi) < tol | mod(timestepping_solution.y(:,7)+2*pi,2*pi) > (2*pi - tol) )        ,1,'first');

time = timestepping_solution.x - timestepping_solution.x(1);
dist = sqrt(sum( (timestepping_solution.y(:,1:3) - repmat(timestepping_solution.y(1,1:3),size(timestepping_solution.y,1),1)).^2 ,2) );
speed = dist ./ time;
time(ind)
speed(ind)
