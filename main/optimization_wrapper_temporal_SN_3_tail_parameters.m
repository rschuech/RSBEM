function [Temporal_SN] = optimization_wrapper_temporal_SN(AR1,AR2,Amp,Lambda,Nlambda, results_file, lock_file_name, results_file_name , sweep_tempfile)

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
    
    Results = temp.Results;
    
    
    if ~isempty(Results)
        already_done = [ [Results.AR1]' [Results.AR2]' [Results.amp]' [Results.lambda]' [Results.nlambda]'  ];
        
        [found, ind] = ismember([AR1 AR2 Amp Lambda Nlambda],already_done,'rows');
        
        if found && ~isempty(Results(ind).taxis) && ~isempty(Results(ind).taxis.temporal.SN)
            Temporal_SN = Results(ind).taxis.temporal.SN;
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
        save(['C:\Users\rudi\Desktop\RD\sweeps\',sweep_tempfile],'new_sweep');
end

settings_inputs_temporal_SN;

[mesh_succeed] = gen_mesh_wrapper(Inputs(1));  %should only be one entry in Inputs, which goes with the tail we need to generate

if ~mesh_succeed
    Temporal_SN = NaN;
    return
end

main;

try
    Temporal_SN = Results(results_ind).taxis.temporal.SN;  %should still be floating around from main.m
catch
    Temporal_SN = NaN;  % probably due to tail intersecting body and most of main.m not running
end
