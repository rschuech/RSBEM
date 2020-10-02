
metric = 'temporal';
% sweep.bodies = [4 0.4; 4 0.5; 5 0.6;];
% sweep.bodies = setdiff(AR1_AR2_data.temporal , AR1_AR2_data.Dm,'rows');  % current temporal is refined, Dm is coarse
% sweep.bodies = [AR1' AR2'];
sweep.bodies = bads;


sb = standardize(sweep.bodies,limits);
params = {'amp','lambda','nlambda'};


for b = 1:size(sweep.bodies,1)
    for p = 1:length(params)
        
        sweep.(params{p})(b,1) = F.(metric).F.(params{p})(sb(b,1),sb(b,2));
    end
end

