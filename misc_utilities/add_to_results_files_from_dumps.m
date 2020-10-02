% add some shat to existing Results by loading corresponding dump files, if
% they exist

% load a Results file here

% dumpfolder = 'C:\Users\rudi\Desktop\RD\efficiency_dumps_pole2pole\';
% dumpfolders_remote = {'X:\efficiency_dumps_pole2pole\' , 'Y:\efficiency_dumps_pole2pole\' , 'Z:\efficiency_dumps_pole2pole\'};

% dumpfolder = 'C:\Users\rudi\Desktop\RD\efficiency_dumps_smooth\';
% dumpfolders_remote = {'X:\efficiency_dumps_smooth\' , 'Y:\efficiency_dumps_smooth\' , 'Z:\efficiency_dumps_smooth\'};

% dumpfolders = {'C:\Users\rudi\Desktop\RD\efficiency_dumps_smooth\',...
%     'X:\efficiency_dumps_smooth\', 'Y:\efficiency_dumps_smooth\' , 'Z:\efficiency_dumps_smooth\',...
%     'C:\Users\rudi\Desktop\RD\efficiency_dumps_refined_optimum\',...
%     'X:\efficiency_dumps_refined_optimum\', 'Y:\efficiency_dumps_refined_optimum\' , 'Z:\efficiency_dumps_refined_optimum\' };

dumpfolders = {'C:\Users\rudi\Desktop\RD\efficiency_dumps_pole2pole\', 'X:\efficiency_dumps_pole2pole\' , 'Y:\efficiency_dumps_pole2pole\' , 'Z:\efficiency_dumps_pole2pole\'};

load_timestepping = false;

phases = linspace(0,2*pi,500);
digits = 16;
workc = 0;

for r = 1:length(Results)
    
    %  curved_rod_AR1_1.1_AR2_0_tail_radius_0.03101752454497_amp_0.4499055576222944_lambda_3.327729624888251_nlambda_1.386798585393612_motorBC_torque_dump.mat
    
    AR1 = Results(r).AR1;  AR2 = Results(r).AR2;  amp = Results(r).amp;  lambda = Results(r).lambda;  nlambda = Results(r).nlambda;
    
    if ~ismember(roundn([AR1 AR2 amp lambda nlambda],-6),roundn([ Best.eff.body  Best.eff.amp Best.eff.lambda Best.eff.nlambda ],-6),'rows')
        Results(r).buckling.force.mean = NaN;  Results(r).buckling.force.max = NaN;
        Results(r).buckling.torque.mean = NaN;  Results(r).buckling.torque.max = NaN;
        Results(r).buckling.torque_ratio.mean = NaN;  Results(r).buckling.torque_ratio.max = NaN;
        continue
    end
    
%     if ~isempty(Results(r).buckling) && ~isnan(Results(r).buckling.torque_ratio.mean)
%         continue
%     end
    
    % if we need something fron freeswim Fourier interp dump
    %     dumpname =   ['curved_rod_AR1_',num2str(AR1),'_AR2_',num2str(AR2),'_tail_radius_0.03101752454497_amp_',num2str(amp,digits),'_lambda_',num2str(lambda,digits),'_nlambda_',num2str(nlambda,digits),'_motorBC_torque_dump.mat'];
    chop_off = 4;
    amp_str = num2str(amp,digits); lambda_str = num2str(lambda,digits);  nlambda_str = num2str(nlambda,digits);
    dumpstr = ['curved_rod_AR1_',num2str(AR1),'_AR2_',num2str(AR2),'_tail_radius_0.03101752454497_amp_',amp_str(1:end-chop_off),'*','_lambda_',lambda_str(1:end-chop_off),'*','_nlambda_',nlambda_str(1:end-chop_off),'*','_motorBC_torque_dump.mat'];
    
    
    
    dump = []; dumpname = [];
    
    for d = 1:length(dumpfolders)
        
        
        temp = dir([dumpfolders{d},dumpstr]);
        if numel(temp) > 1
            error('More than one matching file found');
        end
        if ~isempty(temp)
            dumpname = temp.name;
            
            dump = load([dumpfolders{d}, dumpname]);
            break
            
        end
        
    end
    
    
    if isempty(dump)
        disp('can''t find dump');
        stopa
    end
    
    
    if ~isfield(dump,'Metadata')
        disp('dump doesn''t have Metadata');
        stopa
    end
    
    Metadata = dump.Metadata;
    input = dump.input;
    Mesh = dump.Mesh;
    Solutions = dump.Solutions;
    matrix_props = dump.matrix_props;
    
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
    
    
    
    
    if ~isfield(input.constants,'power')
        input.constants.power = 1E-3;
    end
    
    
    if true || ~isfield(Results,'buckling') || isempty(Results(r).buckling) || isnan(Results(r).buckling.force.mean)
        flagellar_hook_force_torque;  % makes buckling variable
        Results(r).buckling = buckling;
    end
    
    
    
    workc = workc + 1;
    
    doneness =    r/length(Results);
    workfrac = workc / r;
    
    [doneness workfrac]
    
    
    
    
end
