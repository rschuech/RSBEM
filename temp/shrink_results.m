bads = false(1,length(Results));

for r = 1:length(Results)
    if isnan(Results(r).buckling.force.max)
        bads(r) = true;
    end
end