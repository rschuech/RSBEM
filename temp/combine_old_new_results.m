
old_file = 'E:\Hull\Results\current\Results.mat';  %less accurate results to use if we have nothing else
new_file = 'E:\Hull\Results_1.mat';  %new results to prefer over old results



temp = load(old_file);
old_Results = temp.Results;
old_fields = fieldnames(old_Results);

temp = load(new_file);
new_Results = temp.Results;
new_fields = fieldnames(new_Results);

all_Results = new_Results;


old_geoms = [ [old_Results.AR1]'   [old_Results.AR2]'  [old_Results.amp]'  [old_Results.lambda]'  [old_Results.nlambda]'  ];
new_geoms = [ [new_Results.AR1]'   [new_Results.AR2]'  [new_Results.amp]'  [new_Results.lambda]'  [new_Results.nlambda]'  ];

new_bodies = unique(new_geoms(:,1:2),'rows');



for l = 1:length(old_Results)  %each geom result in local Results file
    
    if ~ismember(old_geoms(l,1:2),new_bodies,'rows')  %if this old result is for a body we didn't redo yet
        old_geoms(l,1:2)
        
        if ismember(old_geoms(l,1),[1.125 1.25 1.375 1.625 1.75])
            continue
        end
        
        new_ind = length(all_Results) + 1;
        
        for f = 1:length(old_fields)
            all_Results(new_ind).(old_fields{f}) = old_Results(l).(old_fields{f});
        end
    end
    
end
        
%%
all_geoms = [ [all_Results.AR1]'   [all_Results.AR2]'  [all_Results.amp]'  [all_Results.lambda]'  [all_Results.nlambda]'  ];

all_bodies = unique(all_geoms(:,1:2),'rows');

Results = all_Results;
save('E:\Hull\Results\current\Results_all.mat','Results');