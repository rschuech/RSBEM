function [power_eff] = optimization_wrapper(AR1,AR2,Amp,Lambda,Nlambda, results_file, other_results_files)

fields = {'AR1','AR2','amp','lambda','nlambda','Power_eff'};

if exist(results_file,'file')
    worked = false;
    while ~worked
        try
            temp = load(results_file);
            worked = true;
        catch
            pause(2);
        end
    end
    for f = 1:length(fields)
        for fi = 1:length(temp.Results)
            Results(fi).(fields{f}) = temp.Results(fi).(fields{f});
        end
    end
    for i = 1:length(other_results_files)
        worked = false;
        while ~worked
            try
                temp = load(other_results_files{i});
                worked = true;
            catch
                pause(2);
            end
        end
        clear Results2
        for f = 1:length(fields)
            for fi = 1:length(temp.Results)
                Results2(fi).(fields{f}) = temp.Results(fi).(fields{f});
            end
        end
        Results = [Results Results2];
    end
    
    
    if ~isempty(Results)
        already_done = [ [Results.AR1]' [Results.AR2]' [Results.amp]' [Results.lambda]' [Results.nlambda]'  ];
        
        [found, ind] = ismember([AR1 AR2 Amp Lambda Nlambda],already_done,'rows');
        
        if found
            power_eff = Results(ind).Power_eff;
            return
        end
    end
    
end
% welp, guess we have to run this one

new_sweep.AR1 = AR1;  new_sweep.AR2 = AR2;  new_sweep.amp = Amp;  new_sweep.lambda = Lambda;  new_sweep.nlambda = Nlambda;

switch getenv('computername')
    case 'UBERTOP'
        
        save('E:\Hull\sweeps\new_sweep.mat','new_sweep');
    case {'CFD01','CFD02','CFD03','CFD04'}
        save('C:\Users\rudi\Desktop\RD\sweeps\new_sweep.mat','new_sweep');
end

settings_inputs;

[mesh_succeed] = gen_mesh_wrapper(Inputs(1));  %should only be one entry in Inputs, which goes with the tail we need to generate

if ~mesh_succeed
    power_eff = NaN;
    return
end

main;

power_eff = Results(end).Power_eff;  %should still be floating around from main.m
