

%    POS = [1 7 2 3 5 6]; %Position index of X for repeated values
%    IR  = [1 1 2 2 2 2]; %Corresponding to which index of RV
Results = [];
for rack = 1:4
    temp = load(['E:\Hull\Results\results_CFD0',num2str(rack),'.mat']);
    for i = 1:length(temp.Results)
        temp.Results(i).rack = rack;
    end
    Results = [Results temp.Results];
    
end

%%


geom = [ [Results.AR1]' [Results.AR2]' [Results.amp]' [Results.lambda]' [Results.nlambda]' ];


[RV,NR,POS,IR]=repval(geom);

Repeats = [];
temp = Results(POS(1));

for i= 2:length(POS)
        if IR(i - 1) == IR(i)
            temp(end+1) = Results(POS(i));
        else
            cutoffs = [temp.diff_cutoff];
            [~, best] = min(cutoffs);
            repeats = setdiff(1:length(cutoffs), best);
            for j = 1:length(repeats)
                Repeats(end+1,1).name = temp(repeats(j)).name;
                Repeats(end,1).rack = temp(repeats(j)).rack;
            end
            temp = Results(POS(i));
        end
end
            %%
% save('E:\Hull\Results\repeats.mat','Repeats');
