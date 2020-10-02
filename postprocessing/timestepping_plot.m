

%shat = load([input.paths.dumpfolder,input.paths.fullnamebase,'_timestepping','.mat'],'timestepping_convergence','timestepping_solution');

last_time = shat.timestepping_solution.x(end);

time = linspace(0,last_time,1000);

%xyz = deval(shat.timestepping_solution
x0 = shat.timestepping_solution.y(1:3,1);
dists = sqrt( sum(( shat.timestepping_solution.y(1:3,:) - repmat(x0,1,length(shat.timestepping_solution.x) ) ).^2)) ;

%figure;  plot(shat.timestepping_solution.x, dists, 'o-');

speeds = dists ./ shat.timestepping_solution.x;

figure;  plot(shat.timestepping_solution.x(1:10:end), speeds(1:10:end), 'o-');

