function [power_eff] = optimization_wrapper_efficiency(AR1,AR2,Amp,Lambda,Nlambda ,results_file, lock_file_name, results_file_name , sweep_tempfile, dump_folder)



% obj = @(nlambda) bacterial_tail_arclength([Amp Lambda nlambda]) - max_arclength ;
% % arclength = sqrt(amp^2 + 1/KE^2) * KE*xi    from wolfram helix page, with
% % c = 1/KE
% k = 2*pi./Lambda;  KE = k; % as per Shum et al
% nlambda_guess = max_arclength / sqrt( Amp^2 + 1/KE^2 ) / KE / Lambda;  % using basic helix equation, neglecting variable amp
% Nlambda = fzero(obj, nlambda_guess);





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
        
        if found && ~isempty(Results(ind).Power_eff)
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
        save(['C:\Users\rudi\Desktop\RD\sweeps\',sweep_tempfile],'new_sweep');
end

settings_inputs_medium_forced_freeswim;

[mesh_succeed] = gen_mesh_wrapper(Inputs(1));  %should only be one entry in Inputs, which goes with the tail we need to generate

if ~mesh_succeed
    power_eff = NaN;
    return
end

main;

try
    power_eff = Results(results_ind).Power_eff;  %should still be floating around from main.m
catch
    power_eff = NaN;  % probably due to tail intersecting body and most of main.m not running
end
