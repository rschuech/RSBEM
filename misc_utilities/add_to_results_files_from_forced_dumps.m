% add some shat to existing Results by loading corresponding dump files, if
% they exist

% load a Results file here

dumpfolder = 'C:\Users\rudi\Desktop\RD\body_only_dumps\';
 dumpfolders_remote = {'X:\body_only_dumps\' , 'Y:\body_only_dumps\' , 'Z:\body_only_dumps\'};

load_timestepping = false;

phases = linspace(0,2*pi,1000);
digits = 16;
workc = 0;

for r = 1:length(Results)
    
    %  curved_rod_AR1_1.1_AR2_0_tail_radius_0.03101752454497_amp_0.4499055576222944_lambda_3.327729624888251_nlambda_1.386798585393612_motorBC_torque_dump.mat
    
    AR1 = Results(r).AR1;  AR2 = Results(r).AR2;  amp = Results(r).amp;  lambda = Results(r).lambda;  nlambda = Results(r).nlambda;
    
    
%     dumpname =   ['curved_rod_AR1_',num2str(AR1),'_AR2_',num2str(AR2),'_tail_radius_0.03101752454497_amp_',num2str(amp,digits),'_lambda_',num2str(lambda,digits),'_nlambda_',num2str(nlambda,digits),'_motorBC_torque_dump.mat'];
    dumpname =   ['curved_rod_AR1_',num2str(AR1),'_AR2_',num2str(AR2),'_forced_dump.mat'];
    
    
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
        disp('can''t find dump');
        continue
    end
    
    
    if ~isfield(dump,'Metadata')
        disp('dump doesn''t have Metadata');
        continue
    end
    
    Metadata = dump.Metadata;
    input = dump.input;
    Mesh = dump.Mesh;
%     Solutions = dump.Solutions;
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
    
    
%     if ~isfield(Results,'buckling') || isempty(Results(r).buckling) 
%         flagellar_hook_force_torque;  % makes buckling variable   
%         Results(r).buckling = buckling;
%     end
    
    Results(r).fcoeffs = dump.fcoeffs;
    
    
    [poles] = calc_pole_coords(Metadata);
    poles_vec = poles.head - poles.tail;  poles_vec = poles_vec / sqrt(sum(poles_vec.^2));
    
    % index of rotation principle axis with closest alignment to poles_vec
    % (accounting for both possible directions of each principle axis)
    [~,ind] = min(  min(  acos(  dot( dump.D.rotation.axes , repmat(poles_vec,1,3) , 1))    ,   acos(  dot( -dump.D.rotation.axes , repmat(poles_vec,1,3) , 1))    )   );
    %     clear D
    %     D.rotation.axes = dump.D.rotation.axes(:,[ind setdiff(1:3,ind)]);
    %     D.rotation.diffusivity = dump.D.rotation.diffusivity(:,[ind setdiff(1:3,ind)]);
    
    Results(r).rotational_diffusivity = dump.D.rotation.diffusivity;
      Results(r).rotational_axes = dump.D.rotation.axes;
    %     D.rotation.diffusivity = dump.D.rotation.diffusivity(:,[ind setdiff(1:3,ind)]);
    
    principle_axis = dump.D.rotation.axes(:,ind);
    angle_orig = acos(dot(principle_axis, poles_vec));  angle_reversed = acos(dot(-principle_axis,poles_vec));
    if angle_reversed < angle_orig
        principle_axis = - principle_axis;
    end
    
    
    Results(r).principle_axis_1 =  principle_axis;
    Results(r).principle_axis_1_index = ind;  % we aren't reordering the matrices, just keeping track of the body-oriented axis
    
    workc = workc + 1;
    
    doneness =    r/length(Results);
    workfrac = workc / r;
    
    [doneness workfrac]
    
   
    
    
end
