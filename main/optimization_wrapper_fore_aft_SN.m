function [Fore_Aft_SN] = optimization_wrapper_fore_aft_SN(AR1,AR2,Amp,Lambda,max_arclength ,results_file, lock_file_name, results_file_name , sweep_tempfile)



obj = @(nlambda) bacterial_tail_arclength([Amp Lambda nlambda]) - max_arclength ;
% arclength = sqrt(amp^2 + 1/KE^2) * KE*xi    from wolfram helix page, with
% c = 1/KE
k = 2*pi./Lambda;  KE = k; % as per Shum et al
nlambda_guess = max_arclength / sqrt( Amp^2 + 1/KE^2 ) / KE / Lambda;  % using basic helix equation, neglecting variable amp
Nlambda = fzero(obj, nlambda_guess);





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
        
        if found && ~isempty(Results(ind).taxis) && ~isempty(Results(ind).taxis.fore_aft.SN)
            Fore_Aft_SN = Results(ind).taxis.fore_aft.SN;
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

settings_inputs_fore_aft_SN;

[mesh_succeed] = gen_mesh_wrapper(Inputs(1));  %should only be one entry in Inputs, which goes with the tail we need to generate

if ~mesh_succeed
    Fore_Aft_SN = NaN;
    return
end

main_skip_intersections_check;

try
    Fore_Aft_SN = Results(results_ind).taxis.fore_aft.SN;  %should still be floating around from main.m
catch
    Fore_Aft_SN = NaN;  % probably due to tail intersecting body and most of main.m not running
end
