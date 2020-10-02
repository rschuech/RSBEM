bodies = [1 0; 10 0; 6 0.65; 2 0.5; 10 0.5; 3.5 0.9; 10 0.85;];
clear tails
for b = 1:size(bodies,1)
    
  [~,ind]=  ismember(bodies(b,:),X_Y_unq,'rows');
  
  tails(b,:) = [best_amps(ind) best_lambdas(ind) best_nlambdas(ind)];
  
end


clear new_sweep

new_sweep.AR1 = bodies(:,1)';
new_sweep.AR2 = bodies(:,2)';
new_sweep.amp = tails(:,1)';
new_sweep.lambda = tails(:,2)';
new_sweep.nlambda = tails(:,3)';