


geoms = [ [Results.AR1]' [Results.AR2]' [Results.amp]'  [Results.lambda]' [Results.nlambda]'  ];

bests = [X_Y_unq best_amps best_lambdas best_nlambdas];

[la,lb] = ismember(bests, geoms, 'rows');

Results_best = Results( lb(lb~=0));

geoms_best = [ [Results_best.AR1]' [Results_best.AR2]' [Results_best.amp]'  [Results_best.lambda]' [Results_best.nlambda]'  ];

AR1 = [Results_best.AR1]';
geoms_best = geoms_best( AR1 <= 10, :);
Results_best = Results_best( AR1 <= 10);

Results = Results_best;
save('E:\Hull\Results.mat','Results');

new_sweep.AR1 = geoms_best(:,1)';
new_sweep.AR2 = geoms_best(:,2)';
new_sweep.amp = geoms_best(:,3)';
new_sweep.lambda = geoms_best(:,4)';
new_sweep.nlambda = geoms_best(:,5)';

save('E:\Hull\forced_sweep.mat','new_sweep');

