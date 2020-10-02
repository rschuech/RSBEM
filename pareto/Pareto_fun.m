function [obj] = Pareto_fun(SF, F,limits)

% SF

obj = [  - F.eff.F.metric(standardize(SF,limits)) / F.eff.F.metric(0,0)  ;  - F.temporal.F.metric(standardize(SF,limits)) / F.temporal.F.metric(0,0)  ];