bodies = [1 0; 10 0; 10 0.95; 1.125 0.35; 2 0.635; 3.5 0.945; 5.5 0.985; 5.5 0; 10 0.65; 6 0.65; 3.5 0.25; 5.5 0.15;  10 0.3;  6 0.35;  5.5  0.5;  3 0.5;  10 0.5;];

digits = 16;

for b = 1:size(bodies,1)
    
  [~,ind] =   ismember(bodies(b,:),X_Y_unq,'rows');
  
  tail = [best_amps(ind) best_lambdas(ind) best_nlambdas(ind) ];
  
  dumpnames{b} = ['curved_rod_AR1_',num2str(bodies(b,1)),'_AR2_',num2str(bodies(b,2)),'_tail_radius_0.03101752454497_amp_',num2str(tail(1),digits),'_lambda_',num2str(tail(2),digits),'_nlambda_',num2str(tail(3),digits),'_motorBC_torque_dump.mat'];
  
  
end


save('E:\Hull\dumpnames.mat','dumpnames')

% curved_rod_AR1_1.1_AR2_0.05_tail_radius_0.03101752454497_amp_0.4868462364539689_lambda_3.423686282843095_nlambda_1.305637928757825_motorBC_torque_dump.mat