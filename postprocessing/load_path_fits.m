% folder = 'C:\Users\rudi\Desktop\RD\swept_dumps_archived\';
folder = 'C:\Users\rudi\Desktop\RD\swept_dumps\';
% folder2 = 'E:\Hull\CFD03\swept_dumps\';
if ~exist(results_file)
    save(results_file,'Results');
end


files = dir([folder,'*torque_dump.mat']);
files = {files.name};

% clear Fits Avg_Omega Motor_Torque
%   Fits = struct('T',cell(1,length(files)),'line',cell(1,length(files)),'converged',cell(1,length(files)));
F = length(files);
% Results.file = cell(F,1);

% files = files(randperm(length(files)));

% diff_cutoff = NaN(F,1);
% last_diffs = NaN(F,2);  %last two consecutive percent diffs for avg_speed
% avg_speed = NaN(F,1);
% avg_omega = NaN(F,1);
% geom = NaN(F,5);  %AR1 AR2 amp lambda nlambda
% runtime = NaN(F,1);
% simtime = NaN(F,1);
% motor_torque = NaN(F,1);
% 
% t_converged = NaN(F,1);
% t_goodenough = NaN(F,1);
% cutoff_goodenough = NaN(F,1);
% fucked_diff_cutoff = false(F,1);
% no_file = false(F,1);

speed_tols = [0.1:0.1:5];
Tmax = 2^13;

ppm = ParforProgressStarter2('shat', length(files));

for f = 1:length(files)  %timestepping dumps
    dump_file = files{f};
    %     f/length(files)
    Results(f).filename = files{f};
    
    
    %     Results(i).name = files{i}(1:end-32);  temp =  files{i}; for timestepping files
    Results(f).name = files{f}(1:end-24);  temp =  files{f};
    
    
    Results(f).AR1 = str2double(temp(strfind(temp,'AR1_')+4 : strfind(temp,'_AR2')-1));
    Results(f).AR2 = str2double(temp(strfind(temp,'AR2_')+4 : strfind(temp,'_tail')-1));
    Results(f).amp = str2double(temp(strfind(temp,'amp_')+4 : strfind(temp,'_lambda')-1));
    Results(f).lambda = str2double(temp(strfind(temp,'lambda_')+7 : strfind(temp,'_nlambda')-1));
    Results(f).nlambda = str2double(temp(strfind(temp,'nlambda_')+8 : strfind(temp,'_motorBC')-1));
    
    timestepping_file = [dump_file(1:end-8), 'timestepping.mat'];
    if ~exist([folder, timestepping_file]) %no partial file, load timestepping file if it exists
        stopafro
    end
    temp = load([folder,timestepping_file],'fits','avg_omega','motor_torque');
    temp2 = temp;
    temp2.fits.line.speed = temp2.fits.line.speed(temp2.fits.T <= Tmax);
    temp2.fits.T = temp2.fits.T(temp2.fits.T <= Tmax);
    temp2.fits.converged.diff_cutoff = NaN;  temp2.fits.converged.speed = NaN;
    
    for i = 1:length(speed_tols)
        options = [];  options.mode = 'line';  options.diff_cutoff = speed_tols(i);  options.direction = 'last';
        [temp2.fits] = check_convergence(temp2.fits, options);
        if ~isnan(temp2.fits.converged.speed)
            
            %             t_goodenough(f) = temp2.fits.T( temp2.fits.line.speed == temp2.fits.converged.speed);
            %             cutoff_goodenough(f) = speed_tols(i);
            break
        end
    end
    
    
    Results(f).Avg_Speed = temp2.fits.converged.speed;
    Results(f).diff_cutoff = temp2.fits.converged.diff_cutoff;
    Results(f).t_convergence = temp2.fits.T( temp2.fits.line.speed == temp2.fits.converged.speed);
    if ~isfield(temp,'avg_omega')
        torque_dump = [Results(f).name,'_motorBC_torque_dump.mat'];
        temp4 = load([folder,torque_dump],'interpolant');
        interpolant = temp4.interpolant(end).omega;
        temp.avg_omega = compute_avg_omega(interpolant);
    end
    Results(f).Avg_Omega = temp.avg_omega;
    
    if ~isfield(temp,'motor_torque')
        torque_dump = [Results(f).name,'_motorBC_torque_dump.mat'];
        temp4 = load([folder,torque_dump],'input');
        temp.motor_torque = temp4.input.tail.motor_torque;
    end
    Results(f).Motor_Torque = temp.motor_torque;
    
    
    m = matfile([folder,timestepping_file],'writable',true);
    temp3 = temp.fits;
    temp3.converged = temp2.fits.converged;  %keep orig T and speed values, not truncated version in temp2
    m.fits = temp3;
    m.avg_omega = temp.avg_omega;
    m.motor_torque = temp.motor_torque;
    
    
    ppm.increment(f);
end


delete(ppm);


%%
