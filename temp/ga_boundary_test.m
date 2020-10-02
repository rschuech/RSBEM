min_area = 0.5;

%  [c,ceq] = bnonlcon(inds,AR1AR2,min_area)
% [a , k ] = bobj(inds,AR1AR2)

obj = @(inds) bobj(inds,Pareto_data.observed.species.points);
nonlcon = @(inds) bnonlcon(inds,Pareto_data.observed.species.points,min_area);
nvars = size(Pareto_data.observed.species.points,1);

binds = ga(obj,nvars,[],[],[],[],zeros(1,nvars),ones(1,nvars),nonlcon,1:nvars);