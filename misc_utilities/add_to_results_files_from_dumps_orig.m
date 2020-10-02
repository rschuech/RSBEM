% add some shat to existing Results by loading corresponding dump files, if
% they exist

% load a Results file here

dumpfolder = 'C:\Users\rudi\Desktop\RD\opt_dumps\';
dumpfolders_remote = {'X:\opt_dumps\' , 'Y:\opt_dumps\' , 'Z:\opt_dumps\'}; 

load_timestepping = false;

workc = 0;

for r = 1:length(Results)

    %  curved_rod_AR1_1.1_AR2_0_tail_radius_0.03101752454497_amp_0.4499055576222944_lambda_3.327729624888251_nlambda_1.386798585393612_motorBC_torque_dump.mat
    
    AR1 = Results(r).AR1;  AR2 = Results(r).AR2;  amp = Results(r).amp;  lambda = Results(r).lambda;  nlambda = Results(r).nlambda;
    
    % if we need something fron freeswim Fourier interp dump
    dumpname =   ['curved_rod_AR1_',num2str(AR1),'_AR2_',num2str(AR2),'_tail_radius_0.03101752454497_amp_',num2str(amp,digits),'_lambda_',num2str(lambda,digits),'_nlambda_',num2str(nlambda,digits),'_motorBC_torque_dump.mat'];
    
    dump = [];
    if exist([dumpfolder, dumpname],'file')
         dump = load([dumpfolder, dumpname]);
    else
        for d = 1:length(dumpfolders_remote)
            if exist([dumpfolders_remote{d}, dumpname],'file')
                 dump = load([dumpfolders_remote{d}, dumpname]);
                 break
            end
        end
    end
    if isempty(dump)
        continue
    end
   
    
    if ~isfield(dump,'Metadata')
        continue
    end
    
    Metadata = dump.Metadata;
    input = dump.input;
    
    if load_timestepping
    % if we need something from timestepping dump
    dumpname =   ['curved_rod_AR1_',num2str(AR1),'_AR2_',num2str(AR2),'_tail_radius_0.03101752454497_amp_',num2str(amp,digits),'_lambda_',num2str(lambda,digits),'_nlambda_',num2str(nlambda,digits),'_motorBC_torque_timestepping.mat'];
    
 dump = [];
    if exist([dumpfolder, dumpname],'file')
         dump = load([dumpfolder, dumpname]);
    else
        for d = 1:length(dumpfolders_remote)
            if exist([dumpfolders_remote{d}, dumpname],'file')
                 dump = load([dumpfolders_remote{d}, dumpname]);
                 break
            end
        end
    end
    if isempty(dump)
        continue
    end
    
    
    if ~isfield(dump,'fits')
        continue
    end
    
    fits = dump.fits;
    slope =   fits.line.slope(:,fits.line.speed == fits.converged.speed);
    
    
    end
    
    
    Results(r).path_slope = slope;
    
    Results(r).pole_separation = calc_pole_separation(input, Metadata, fits);
    
        if ~isfield(input.constants,'power')
            input.constants.power = 1E-3;
        end
        
         Results(r).Adj_Speed = sqrt( Results(r).Avg_Speed.^2 .* input.constants.power  ./ Results(r).Avg_Power);
        
           Results(r).Dm = 1/3 * Results(r).Adj_Speed.^2 * Results(r).tau_a;
  swimming_axis = Results(r).path_slope ./ sqrt(sum(Results(r).path_slope.^2));   %3 x 1 each between -1, 1
     Results(r).error_angle = acos(dot( Results(r).principle_axis_1, swimming_axis)) * 180/pi;  %how big is erroneous angle between rotational diffusion principle axis and avg swimming axis
  Results(r).taxis.temporal.SN = Results(r).Adj_Speed * sqrt(Results(r).tau_a);
        Results(r).taxis.temporal.ability = Results(r).taxis.temporal.SN * Results(r).Adj_Speed;
  Results(r).taxis.fore_aft.SN = Results(r).pole_separation *  sqrt(Results(r).tau_a);
        Results(r).taxis.fore_aft.ability =  Results(r).taxis.fore_aft.SN * Results(r).Adj_Speed;
  
        
        
        
    
    workc = workc + 1;
    
    
     doneness =    r/length(Results);
     workfrac = workc / r;
     
     [doneness workfrac]
    
    
    
    
end

stopa
%%
temp = [Inputs.body];  temp2 = [Inputs.tail];
inputs_geom = [  [temp.AR]'  [temp2.amp]' [temp2.lambda]' [temp2.nlambda]'  ];

for results_ind = 1:length(Results)
    
    results_geom = [Results(results_ind).AR1  Results(results_ind).AR2  Results(results_ind).amp  Results(results_ind).lambda  Results(results_ind).nlambda];
    
    [~,lib] = ismember(results_geom, inputs_geom,'rows');
    input = Inputs(lib);
    
         Results(results_ind).Adj_Speed = sqrt( Results(results_ind).Avg_Speed.^2 .* input.constants.power  ./ Results(results_ind).Avg_Power);
        
        Results(results_ind).pole_separation = calc_pole_separation(input, Metadata, fits);
        
           Results(results_ind).Dm = 1/3 * Results(results_ind).Adj_Speed.^2 * Results(results_ind).tau_a;
  swimming_axis = Results(results_ind).path_slope ./ sqrt(sum(Results(results_ind).path_slope.^2));   %3 x 1 each between -1, 1
     Results(results_ind).error_angle = acos(dot( Results(results_ind).principle_axis_1, swimming_axis)) * 180/pi;  %how big is erroneous angle between rotational diffusion principle axis and avg swimming axis
  Results(results_ind).taxis.temporal.SN = Results(results_ind).Adj_Speed * sqrt(Results(results_ind).tau_a);
        Results(results_ind).taxis.temporal.ability = Results(results_ind).taxis.temporal.SN * Results(results_ind).Adj_Speed;
  Results(results_ind).taxis.fore_aft.SN = Results(results_ind).pole_separation *  sqrt(Results(results_ind).tau_a);
        Results(results_ind).taxis.fore_aft.ability =  Results(results_ind).taxis.fore_aft.SN * Results(results_ind).Adj_Speed;
  
end