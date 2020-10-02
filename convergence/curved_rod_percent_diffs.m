%% Fourier interp max rel diff
reldiffs = [2 1 0.5 0.25 0.05 0.01];  % Fourier interp tolerance

clear power_effs

for i = 1:length(reldiffs)

load(['C:\Users\rudi\Desktop\RD\Results\Results_convergence_max_rel_diff_',num2str(reldiffs(i)),'.mat']);

power_effs(i,:) = [Results.Power_eff];

end

percent_diffs = abs( ( power_effs - repmat(power_effs(end,:), size(power_effs,1),1) ) )  ./ repmat(power_effs(end,:), size(power_effs,1),1)  * 100


% 7 bodies tested were [1 0; 10 0; 6 0.65; 2 0.5; 10 0.5; 3.5 0.9; 10 0.85;]
% with corresponding best tails

% max rel diff = 0.5 gives all % diffs < 0.025%
% seems to usually require 8 phase evals


%% adaptive integration tols
clear power_effs tau

load('C:\Users\rudi\Desktop\RD\Results\Results_convergence_max_rel_diff_0.5_reltols_1.mat');
power_effs(1,:) = [Results.Power_eff];
tau(1,:) = [Results.tau_a];

load('C:\Users\rudi\Desktop\RD\Results\Results_convergence_max_rel_diff_0.5_reltols_0.25.mat');
power_effs(2,:) = [Results.Power_eff];
tau(2,:) = [Results.tau_a];

percent_diffs = abs( ( power_effs - repmat(power_effs(end,:), size(power_effs,1),1) ) )  ./ repmat(power_effs(end,:), size(power_effs,1),1)  * 100
percent_diffs = abs( ( tau - repmat(tau(end,:), size(tau,1),1) ) )  ./ repmat(tau(end,:), size(tau,1),1)  * 100

% all % diffs are < 0.03%

