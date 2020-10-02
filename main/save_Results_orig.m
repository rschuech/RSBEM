
clear fid
worked = false;
tic
while ~worked
    
    if exist( [input.paths.global_lock_file],'file' )
        disp(['global lock file found, paused for ',num2str(toc/60),' min']);
        pause(20);
        continue
    end
    
    worked = true;
    fid = fopen( [input.paths.results_folder,input.paths.lock_file] , 'w' );  fclose(fid);
    
end



worked = false;
tic
while ~worked
    try
        temp = load([input.paths.results_folder,input.paths.results_file]);
        Results = temp.Results;
        worked = true;
    catch
        disp(['can''t load local lock file, paused for ',num2str(toc/60),' min']);
        pause(5)
    end
end
     


if ~isempty(Results)
    geoms = [ [Results.AR1]'   [Results.AR2]'  [Results.amp]'  [Results.lambda]'  [Results.nlambda]'  ];
    geom_current = [   input.body.AR(1)   input.body.AR(2)  input.tail.amp  input.tail.lambda  input.tail.nlambda  ];
    
    [~,results_ind] =      ismember(geom_current, geoms, 'rows');
else
    results_ind = 0;
end

if results_ind == 0 % don't have this geom in Results yet, add a new entry
    results_ind = length(Results) + 1;
    % else, ind stays as above
end

Results(results_ind).name = [input.paths.namebase.body, '_', input.paths.namebase.tail];
Results(results_ind).AR1 = input.body.AR(1);
Results(results_ind).AR2 = input.body.AR(2);
Results(results_ind).amp = input.tail.amp;
Results(results_ind).lambda = input.tail.lambda;
Results(results_ind).nlambda = input.tail.nlambda;


if       strcmp(input.problemtype,'freeswim')  && input.do_timestepping
    
    Results(results_ind).Avg_Speed = fits.converged.speed;
    Results(results_ind).path_slope = fits.line.slope(:,fits.line.speed == fits.converged.speed);
    Results(results_ind).diff_cutoff = fits.converged.diff_cutoff;
    Results(results_ind).t_convergence = fits.T( fits.line.speed == fits.converged.speed);
    Results(results_ind).Avg_Omega = avg_omega;
    Results(results_ind).Motor_Torque = input.tail.motor_torque;
    
    Results(results_ind).Avg_Power = Results(results_ind).Avg_Omega * Results(results_ind).Motor_Torque;  %one of these should be constant and the other varying unless I eventually implement a constant power motor condition
    
    Results(results_ind).Torque_eff = 8*pi*input.constants.mu*input.body.sphererad^2*Results(results_ind).Avg_Speed / Results(results_ind).Motor_Torque;
    
    Results(results_ind).Power_eff = 6*pi*input.constants.mu*input.body.sphererad*Results(results_ind).Avg_Speed^2 / Results(results_ind).Avg_Power;
    
    Results(results_ind).Adj_Speed = sqrt( Results(results_ind).Avg_Speed.^2 .* input.constants.power  ./ Results(results_ind).Avg_Power);
    
    Results(results_ind).pole_separation = calc_pole_separation(input, Metadata, fits);
    
end





if       strcmp(input.problemtype,'forced')
    Results(results_ind).tau_a = 1/ (D.rotation.diffusivity(2,2) +  D.rotation.diffusivity(3,3));  %timescale for loss of orientation around first principle axis, which is hopefully close to swimming axis
    
    
    rotation_axis = D.rotation.axes(:,1) / sqrt(sum(D.rotation.axes(:,1).^2));
    if rotation_axis(1) < 0 %principle axis happens to point backward
        rotation_axis = - rotation_axis;
    end
    
    Results(results_ind).principle_axis_1 = rotation_axis;  %the principle axis of rotational diffusion closest to the swimming direction
    
end




% needs both freeswim and forced
% outputs - update all these after
% every run
if isfield(Results, 'Adj_Speed') && isfield(Results,'tau_a')
    Results(results_ind).Dm = 1/3 * Results(results_ind).Adj_Speed.^2 * Results(results_ind).tau_a;
else
    Results(results_ind).Dm = [];
end

if isfield(Results,'path_slope')
    swimming_axis = Results(results_ind).path_slope ./ sqrt(sum(Results(results_ind).path_slope.^2));   %3 x 1 each between -1, 1
else
    swimming_axis = [];
end
if isfield(Results,'principle_axis_1') && ~isempty(Results(results_ind).principle_axis_1) && ~isempty(swimming_axis)
    Results(results_ind).error_angle = acos(dot( Results(results_ind).principle_axis_1, swimming_axis)) * 180/pi;  %how big is erroneous angle between rotational diffusion principle axis and avg swimming axis
else
    Results(results_ind).error_angle = [];
end

if isfield(Results,'Adj_Speed') && isfield(Results,'tau_a')
    Results(results_ind).taxis.temporal.SN = Results(results_ind).Adj_Speed * sqrt(Results(results_ind).tau_a);
    Results(results_ind).taxis.temporal.ability = Results(results_ind).taxis.temporal.SN * Results(results_ind).Adj_Speed;
else
    Results(results_ind).taxis.temporal.SN =  [];
    Results(results_ind).taxis.temporal.ability = [];
end

if isfield(Results,'pole_separation') && isfield(Results,'tau_a') &&  isfield(Results,'Adj_Speed')
    Results(results_ind).taxis.fore_aft.SN = Results(results_ind).pole_separation *  sqrt(Results(results_ind).tau_a);
    Results(results_ind).taxis.fore_aft.ability =  Results(results_ind).taxis.fore_aft.SN * Results(results_ind).Adj_Speed;
else
    Results(results_ind).taxis.fore_aft.SN = [];
    Results(results_ind).taxis.fore_aft.ability =  [];
end


tic
while true
    try
        save([input.paths.results_folder,input.paths.results_file] , 'Results'  );
        break
    catch
        disp(['can''t save local results file, paused for ',num2str(toc/60),' min']);
        pause(5)
    end
end




delete([input.paths.results_folder,input.paths.lock_file]);