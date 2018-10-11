
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
        disp(['can''t load local Results file, paused for ',num2str(toc/60),' min']);
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
    Results(results_ind).avg_swimming_axis =  fits.avg_swimming_axis;
    Results(results_ind).diff_cutoff = fits.converged.diff_cutoff;
    Results(results_ind).t_convergence = fits.T( fits.line.speed == fits.converged.speed);
    Results(results_ind).Avg_Omega = avg_omega;
    Results(results_ind).Motor_Torque = input.tail.motor_torque;
    
    Results(results_ind).Avg_Power = Results(results_ind).Avg_Omega * Results(results_ind).Motor_Torque;  %one of these should be constant and the other varying unless I eventually implement a constant power motor condition
    
    Results(results_ind).Torque_eff = 8*pi*input.constants.mu*input.body.sphererad^2*Results(results_ind).Avg_Speed / Results(results_ind).Motor_Torque;
    
    Results(results_ind).Power_eff = 6*pi*input.constants.mu*input.body.sphererad*Results(results_ind).Avg_Speed^2 / Results(results_ind).Avg_Power;
    
    Results(results_ind).Adj_Speed = sqrt( Results(results_ind).Avg_Speed.^2 .* input.constants.power  ./ Results(results_ind).Avg_Power);
    
    Results(results_ind).pole_separation = calc_pole_separation(input, Metadata, fits.avg_swimming_axis);
    
end





if       strcmp(input.problemtype,'forced')
    
    Results(results_ind).rotational_diffusivity = D.rotation.diffusivity;
    
    
    rotation_axis = D.rotation.axes(:,1) / sqrt(sum(D.rotation.axes(:,1).^2));
    if rotation_axis(1) < 0 %principle axis happens to point backward
        rotation_axis = - rotation_axis;
    end
    
    if ~input.include_tail && isequal( input.body.AR' , [1 0] ) %in degenerate case of an exact sphere with no tail, code computes random principle axis, need to override
        Results(results_ind).principle_axis_1 = [1 0 0]';
    else
        
        Results(results_ind).principle_axis_1 = rotation_axis;  %the principle axis of rotational diffusion closest to the swimming direction
    end
    
    
end




% needs both freeswim and forced
% outputs - update all these after
% every run


if input.swimming_axis_override  % override avg_swimming_axis with those from optimized Temporal SN?
    % for Fore-Aft S/N runs only, hardcode avg_swimming_axis as same as optimized Temporal SN since, for straight tails, actual
    % "swimming path" is some random direction at extremely slow speed
    load(input.swimming_axis_override_file);
    ind = ismember(input.body.AR',avg_swimming_axes.AR1_AR2,'rows');
    
    Results(results_ind).avg_swimming_axis = avg_swimming_axes.avg_swimming_axis(ind,:)';
    
    Results(results_ind).pole_separation = calc_pole_separation(input, Metadata,  Results(results_ind).avg_swimming_axis);
    
end


if isfield(Results,'avg_swimming_axis')
    swimming_axis = Results(results_ind).avg_swimming_axis ./ sqrt(sum(Results(results_ind).avg_swimming_axis.^2));   %3 x 1 each between -1, 1
else
    swimming_axis = [];
end

if isfield(Results,'principle_axis_1') && ~isempty(Results(results_ind).principle_axis_1) && ~isempty(swimming_axis) && isfield(Results,'rotational_diffusivity') && ~isempty(Results(results_ind).rotational_diffusivity)
    Results(results_ind).error_angle = acos(dot( Results(results_ind).principle_axis_1, swimming_axis)) * 180/pi;  %how big is erroneous angle between rotational diffusion principle axis and avg swimming axis
   
    
 principle_axis_1 = Results(results_ind).principle_axis_1;
    rotvec = crossprod( principle_axis_1 , swimming_axis);  rotvec = rotvec / sqrt(sum(rotvec.^2));
    angle = acos( dot(swimming_axis, principle_axis_1) );
    rotmat = rotate_arbitrary_vector( rotvec, angle);
    D_aligned = rotmat * Results(results_ind).rotational_diffusivity * rotmat';  % that's how you rotate a tensor
    
    
    
 fun = @(angle) - rotate_off_axis_tau(angle, swimming_axis, D_aligned);

guesses = linspace(0,2*pi,10);  angle = NaN(size(guesses));  max_tau = angle;  % 10 evenly spaced initial guesses seems quite enough (3 would probably do)
for nn = 1:length(guesses)
[angle(nn),temp] = fminsearch(fun,guesses(nn));
max_tau(nn) = -temp;
end
[max_tau,ind] = max(max_tau);
angle = angle(ind);



% changed to above to match method for body only
%     fun = @(angle) rotate_off_axis_tau(angle, swimming_axis, D_aligned);
%     limits = [0 2*pi];
%     
%     tau_mean = 1/diff(limits) * integral(fun, limits(1),limits(2),'reltol',1E-9,'abstol',1E-12);
    
    
    Results(results_ind).tau_a = max_tau;  %timescale for loss of orientation around swimming axis, ignoring fact that this is not the principle axis and there are cross terms.
  
    
else
    Results(results_ind).error_angle = [];
    Results(results_ind).tau_a = [];
end

if isfield(Results, 'Adj_Speed') && isfield(Results,'tau_a')
    Results(results_ind).Dm = 1/3 * Results(results_ind).Adj_Speed.^2 * Results(results_ind).tau_a;
else
    Results(results_ind).Dm = [];
end


if isfield(Results,'Adj_Speed') && isfield(Results,'tau_a')
    Results(results_ind).taxis.temporal.SN = Results(results_ind).Adj_Speed * (Results(results_ind).tau_a)^(3/2);
%     Results(results_ind).taxis.temporal.ability = Results(results_ind).taxis.temporal.SN * Results(results_ind).Adj_Speed;
else
    Results(results_ind).taxis.temporal.SN =  [];
%     Results(results_ind).taxis.temporal.ability = [];
end

if isfield(Results,'pole_separation') && isfield(Results,'tau_a') 
    Results(results_ind).taxis.fore_aft.SN = Results(results_ind).pole_separation *  sqrt(Results(results_ind).tau_a);
%     Results(results_ind).taxis.fore_aft.ability =  Results(results_ind).taxis.fore_aft.SN * Results(results_ind).Adj_Speed;
else
    Results(results_ind).taxis.fore_aft.SN = [];
%     Results(results_ind).taxis.fore_aft.ability =  [];
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