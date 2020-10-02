
Results_shrunk = [];


bodies = unique([ [Results.AR1]' [Results.AR2]'  ],'rows');

bods_effs = [ [Results.AR1]' [Results.AR2]'  [Results.Power_eff]'];

for b = 1:size(bodies,1)
    body = bodies(b,:);
    
    inds = find(ismember(bods_effs(:,1:2),body,'rows'));
    
    effs = bods_effs(inds,3); %all effs for this body
    [~,ind] = max(effs); %single best eff for this body
    results_ind = inds(ind);  %where in Results is the max eff for this body
    
    if b == 1  %shatlab doesn't like to index into the first one
        Results_shrunk = Results(results_ind);
       
    else
        Results_shrunk(b) = Results(results_ind);
    end
end