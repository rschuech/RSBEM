

% objfun = @(SF) [  - F.eff.F.metric(standardize(SF,limits)) / F.eff.F.metric(0,0)  ;  - F.temporal.F.metric(standardize(SF,limits)) / F.temporal.F.metric(0,0)  ];
% 
% options = optimoptions('gamultiobj','FunctionTolerance',1E-6,'ParetoFraction',0.2,'HybridFcn',@fgoalattain,'PopulationSize',1E16,'UseParallel',true);
% sol = gamultiobj(objfun, 2, [],[],[],[],[1 0],[10 0.8]);
% 
% % sol = unstandardize(sol, limits);
% 
% figure(392);
% plot(sol(:,1),sol(:,2),'o'); grid on;  xlim([1 10]); ylim([0 1]);

%%

objfun = @(SF) - F.eff.F.metric(standardize(SF,limits)) / F.eff.F.metric(0,0);
options = optimset('display','none','tolx',1E-6,'tolfun',1E-6);
[optima.eff.SF, optima.eff.metric] = fminsearch(objfun, [1 , 0],options);

objfun = @(SF) - F.temporal.F.metric(standardize(SF,limits)) / F.temporal.F.metric(0,0);
options = optimset('display','none','tolx',1E-6,'tolfun',1E-6);
[optima.temporal.SF, optima.temporal.metric] = fminsearch(objfun, [9 , 0],options);

% objfun = @(SF) [  - F.eff.F.metric(standardize(SF,limits)) / F.eff.F.metric(0,0)  ;  - F.temporal.F.metric(standardize(SF,limits)) / F.temporal.F.metric(0,0)  ];

objfun = @(SF) Pareto_fun(SF, F, limits);

goal = [optima.eff.metric , optima.temporal.metric];
options = optimoptions('fgoalattain','MaxFunctionEvaluations',1E6,'Display','final-detailed','MaxIterations',1E6);
N = 10;
weights = [ linspace(0,1,N)' (1-linspace(0,1,N))'];
guesses = [ linspace(10,1,N)'  repmat(0.5,N,1)  ];
lb = [4 0.0];  ub = [10 0.9];
clear sol fval
for w = 1:size(weights,1)
[sol(w,:) , fval(w,:)] = fgoalattain(objfun, guesses(w,:) , goal , weights(w,:) ,[],[],[],[],[lb],[ub],[],options);
end


