function [power_eff] = optimization_wrapper_redo_timestepping(AR1,AR2,Amp,Lambda,Nlambda, results_file, lock_file_name, results_file_name , sweep_tempfile)

digits = 16;
% old simulations with broken timestepping but working interp
dump_dirs = {'C:\Users\rudi\Desktop\RD\interped_guesses_dumps\'  ,  'X:\interped_guesses_dumps\' , 'Y:\interped_guesses_dumps\' , 'Z:\interped_guesses_dumps\'};

dump_name = ['curved_rod_AR1_',num2str(AR1),'_AR2_',num2str(AR2),'_tail_radius_0.03101752454497_amp_',num2str(Amp,digits),'_lambda_',num2str(Lambda,digits),...
    '_nlambda_',num2str(Nlambda,digits),'_motorBC_torque_dump.mat'];


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
    
    Results = temp.Results;
    
    
    if ~isempty(Results)
        already_done = [ [Results.AR1]' [Results.AR2]' [Results.amp]' [Results.lambda]' [Results.nlambda]'  ];
        
        [found, ind] = ismember([AR1 AR2 Amp Lambda Nlambda],already_done,'rows');
        
        if found && ~isempty(Results(ind).Power_eff) % we did this run before
            
            
            
            power_eff = Results(ind).Power_eff;
            return
        end
    end
end

% wasn't in current Results, now check broken Results, maybe can just redo
% timestepping
% already_done = [ [Results_broken.AR1]' [Results_broken.AR2]' [Results_broken.amp]' [Results_broken.lambda]' [Results_broken.nlambda]'  ];

%         [found, ind] = ismember([AR1 AR2 Amp Lambda Nlambda],already_done,'rows');

%         if found && ~isempty(Results(ind).Power_eff) % we did this run before, Power_eff is there but it's wrong
% try to find the corresponding dump file
for d = 1:length(dump_dirs)
    if exist([dump_dirs{d},dump_name],'file')
        load([dump_dirs{d},dump_name],'y0','best_interpolant','interpolant','input','Mesh','Metadata');
        
        % update certain inputs
        input.paths.dumpfolder = dump_dirs{d}; %overwrite broken timestepping dump in original location
        input.paths.results_file = results_file_name;
        input.accuracy.interpolation.max_rel_diff = 0.5 ;
        input.accuracy.interpolation.vector_normalization = true;
        input.accuracy.interpolation.max_iter = 5;
        input.accuracy.timestepping.T_interrogate = 2.^(-2:9);
        input.accuracy.timestepping.diff_tols = [0.05 0.1:0.1:5];
        
        
        if ~exist([input.paths.results_folder,input.paths.results_file],'file')
            Results = [];
            save([input.paths.results_folder,input.paths.results_file],'Results');
        end
        
        
        
        dump.y0 = y0;  dump.best_interpolant = best_interpolant;
        
        timestepping_solution0 = [];  refpoint0 = Mesh(1).refpoints(:,1);
        [fits, timestepping_solution, timings] = timestepping(input,dump,timestepping_solution0,refpoint0);  %timestepping dump file is saved internally, contains timestepping_solution and fits
        
        
        interpolant = interpolant(end).omega;  %last and most accurate interpolant for tail rotation rate omega
        
        [avg_omega] = compute_avg_omega(interpolant);
        
        %in addition to storing in memory for later inclusion into aggregate results file, immediately save avg_omega and motor_torque inside timestepping dump file
        m = matfile([input.paths.dumpfolder,input.paths.namebase.full,'_timestepping','.mat'],'Writable',true);
        m.avg_omega = avg_omega;
        m.motor_torque = input.tail.motor_torque;
        
        fits.avg_swimming_axis = calc_avg_swimming_axis(fits, timestepping_solution, Mesh); %outputs direction of avg swimming direction in body frame
        
        m.fits = fits;
        
        
        save_Results;
        
        power_eff = Results(results_ind).Power_eff;
        return
    end
end
disp('old dump file not found, starting from scratch')



% welp, guess we have to run this one

new_sweep.AR1 = AR1;  new_sweep.AR2 = AR2;  new_sweep.amp = Amp;  new_sweep.lambda = Lambda;  new_sweep.nlambda = Nlambda;

switch getenv('computername')
    case 'UBERTOP'
        
        save('E:\Hull\sweeps\new_sweep.mat','new_sweep');
    case {'CFD01','CFD02','CFD03','CFD04'}
        save(['C:\Users\rudi\Desktop\RD\sweeps\',sweep_tempfile],'new_sweep');
end

settings_inputs;

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
