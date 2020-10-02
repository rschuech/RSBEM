
goals = {'Construction_Ease','Power_eff','Dm','temporal_ability','fore_aft_ability'};

ranks = [];
for g  = 1:length(goals)
    ranks = cat(3,ranks,Pareto_data.(goals{g}).rank);
end

optimal = true(size(Pareto_data.Power_eff.rank));

for i = 1:numel(Pareto_data.Power_eff.rank)
    i
    [r1,c1] = ind2sub(size(Pareto_data.Power_eff.rank),i);
    for ii = 1:numel(Pareto_data.Power_eff.rank)
        
        if ~optimal(ii)
            continue
        end
        [r2,c2] = ind2sub(size(Pareto_data.Power_eff.rank),i);
        % for g = 1:length(goals)
        if all( ranks(r,c,:)  < (ranks(r2,c2,:)) )  %another point exists that's better in all goals than this point
            optimal(i) = false;
            
            break;
        end
    end
    sum(~optimal(:))
end