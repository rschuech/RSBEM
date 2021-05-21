yprime = @(t,y) 1 + t.*sin(t.*y);

options = odeset('Refine',[1],'OutputFcn',@output_fun,'MaxStep',0.2,'RelTol',Inf,'AbsTol',Inf);

sol = ode45(yprime, [1 5], 0,options);

diff(sol.x)

figure(345);  plot(sol.x,sol.y,'o-')
figure(342); plot(diff(sol.x),'o');


function status = output_fun(t,y,flag)
t;
status = 0;
end