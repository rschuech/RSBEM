



median_pt = median( vertcat( species_data.mean_unweighted ) );
median_pt_std = standardize(median_pt,limits);

eff_tail = [F.eff.F.amp(median_pt_std) F.eff.F.lambda(median_pt_std)  F.eff.F.nlambda(median_pt_std) ];

Dm_tail = [F.Dm.F.amp(median_pt_std) F.Dm.F.lambda(median_pt_std)  F.Dm.F.nlambda(median_pt_std) ];

temporal_tail = [F.temporal.F.amp(median_pt_std) F.temporal.F.lambda(median_pt_std)  F.temporal.F.nlambda(median_pt_std) ];

clear sweep
sweep.bodies = Best.eff.body;
sweep.amp = repmat(eff_tail(1),size(sweep.bodies,1),1);
sweep.lambda = repmat(eff_tail(2),size(sweep.bodies,1),1);
sweep.nlambda = repmat(eff_tail(3),size(sweep.bodies,1),1);
save('C:\Hull\sweeps\median_obs_body_eff_opt_tail_sweep.mat','sweep');

clear sweep
sweep.bodies = Best.Dm.body;
sweep.amp = repmat(Dm_tail(1),size(sweep.bodies,1),1);
sweep.lambda = repmat(Dm_tail(2),size(sweep.bodies,1),1);
sweep.nlambda = repmat(Dm_tail(3),size(sweep.bodies,1),1);
save('C:\Hull\sweeps\median_obs_body_Dm_opt_tail_sweep.mat','sweep');

clear sweep
sweep.bodies = Best.temporal.body;
sweep.amp = repmat(temporal_tail(1),size(sweep.bodies,1),1);
sweep.lambda = repmat(temporal_tail(2),size(sweep.bodies,1),1);
sweep.nlambda = repmat(temporal_tail(3),size(sweep.bodies,1),1);
save('C:\Hull\sweeps\median_obs_body_temporal_opt_tail_sweep.mat','sweep');



