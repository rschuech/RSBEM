% folder = 'C:\Users\rudi\Desktop\RD\swept_dumps_archived\';
folder = 'C:\Users\rudi\Desktop\RD\swept_dumps\';
% folder2 = 'E:\Hull\CFD03\swept_dumps\';
results_file = 'C:\Users\rudi\Desktop\RD\results_CFD04.mat';
Results = [];
if ~exist(results_file)
    save(results_file,'Results');
end


files = dir([folder,'*torque_dump.mat']);
files = {files.name};

% clear Fits Avg_Omega Motor_Torque
%   Fits = struct('T',cell(1,length(files)),'line',cell(1,length(files)),'converged',cell(1,length(files)));
F = length(files);
% Results.file = cell(F,1);

files = files(randperm(length(files)));

diff_cutoff = NaN(F,1);
last_diffs = NaN(F,2);  %last two consecutive percent diffs for avg_speed
avg_speed = NaN(F,1);
avg_omega = NaN(F,1);
geom = NaN(F,5);  %AR1 AR2 amp lambda nlambda
runtime = NaN(F,1);
simtime = NaN(F,1);
motor_torque = NaN(F,1);


ppm = ParforProgressStarter2('shat', length(files));

parfor f = 1:length(files)  %timestepping dumps
    dump_file = files{f};
    %     f/length(files)
    
    %   if ~ismember(f,B)
    %         continuehhh
    %     end
    %     if isstruct(Fits(f).converged)
    %         continue
    %     end
    %     Results.file(f) = dump_file;
    
    
    
    %     try
    
    %  interp_file = [ timestepping_file(1:end-16) 'dump.mat'];
    
    %         if ~exist([folder,interp_file])
    %             continue
    %         end
    
    %         newfile = [timestepping_file(1:end-16) 'helix_fit.mat'];
    %
    %     if exist([folder2,newfile],'file')
    %         continue
    %     end
    
    timestepping_file = [dump_file(1:end-8), 'timestepping.mat'];
    if exist([folder, timestepping_file]) %no partial file, load timestepping file if it exists
        temp2 = load([folder,timestepping_file]);
        if temp2.fits.converged.diff_cutoff == 0.1 || temp2.timestepping_solution.x(end) >= 2^13 * 0.99
            ppm.increment(f);
            continue
        end
    end
    
    temp = load([folder,dump_file]);
    
    input = temp.input;
    
    
    
    if isfield(input.paths,'namebase')
        name = input.paths.namebase.full;
    else
        name = input.paths.fullnamebase;
    end
    
    
    %     if exist( [input.paths.dumpfolder,name,'_timestepping','.mat'] )
    %           ppm.increment(f);
    %         continue
    %     end
    
    
    
    input.accuracy.timestepping.initialstep = 1E-6;
    input.accuracy.timestepping.reltol = 1E-7;
    input.accuracy.timestepping.abstol = 1E-14;
    input.accuracy.timestepping.T_interrogate = 2.^[-2:13];
    
    dump = [];
    dump.y0 = temp.y0;  dump.best_interpolant = temp.best_interpolant;
    
    
    %first look for partial timestepping file, which is likely further
    %along than orig timestepping dump
    
%     partial_file = [name,'_partial.mat'];
%     if exist([folder,partial_file])
%         temp2 = load([folder,partial_file]);  %timestepping_solution and fits
%         timestepping_solution0 = temp2.timestepping_solution;
%     elseif exist([folder, name, '_timestepping.mat']) %no partial file, load timestepping file if it exists
%         temp2 = load([folder,name,'_timestepping.mat']);
%         timestepping_solution0 = temp2.timestepping_solution;
%     else
%         timestepping_solution0 = [];
%     end
      timestepping_solution0 = [];
    [fits, timestepping_solution, timings] = timestepping(input,dump,timestepping_solution0);  %timestepping dump file is saved internally, contains timestepping_solution and fits
    
    diff_cutoff(f) = fits.converged.diff_cutoff;
    last_diffs(f,:) = abs(  diff(fits.line.speed(end-2:end)) ./ fits.line.speed(end-2:end-1) ) * 100;  %last two consecutive percent diffs for avg_speed
    avg_speed(f) = fits.converged.speed;
    
    % get avg omega from interpolant in interpolant dump file
    %         temp = load([folder,interp_file],'interpolant','input');
    interpolant = temp.interpolant(end).omega;  %last and most accurate interpolant for tail rotation rate omega
    
    [avg_omega(f)] = compute_avg_omega(interpolant);
    
    %in addition to storing in memory for later inclusion into aggregate results file, immediately save avg_omega and motor_torque inside timestepping dump file
    m = matfile([input.paths.dumpfolder,name,'_timestepping','.mat'],'Writable',true);
    m.avg_omega = avg_omega(f);
    m.motor_torque = input.tail.motor_torque;
    
    geom(f,:) = [ temp.input.body.AR' temp.input.tail.amp temp.input.tail.lambda temp.input.tail.nlambda];  %AR1 AR2 amp lambda nlambda
    
    runtime(f) = timings.timestepping;
    simtime(f) = timestepping_solution.x(end);
    motor_torque(f) = input.tail.motor_torque;
    
    
    
    
    %     try
    %     temp = load(results_file);
    %
    %
    %     Results = temp.Results;
    
    %     clear results
    %     results.fit = Fits(i);
    %     results.Avg_Omega = Avg_Omega(i);
    %     results.Motor_Torque = Motor_Torque(i);
    %
    %         Avg_Omega(f) = avg_omega;
    %
    %         Motor_Torque(f) = temp.input.tail.motor_torque;
    % save([folder2,newfile],'fits');
    
    
    
    %     catch
    %         Fits(f).T = NaN;
    %         Fits(f).line.intercept = NaN;
    %         Fits(f).line.slope = NaN;
    %         Fits(f).line.speed = NaN;
    %         Fits(f).line.iters = NaN;
    %         Fits(f).line.flag = NaN;
    %
    %         Fits(f).converged.diff_cutoff = NaN;
    %         Fits(f).converged.speed = NaN;
    %
    %         Avg_Omega(f) = NaN;
    %         Motor_Torque(f) = NaN;
    %
    %     end
    
    ppm.increment(f);
    
    
    
end

delete(ppm);

stopafra

%%

clear Results
for i = 1:length(Fits)
    Results(i).fit = Fits(i);
    Results(i).Avg_Omega = Avg_Omega(i);
    Results(i).Motor_Torque = Motor_Torque(i);
    Results(i).name = files{i}(1:end-32);  temp =  files{i};
    Results(i).filename = files{i};
    Results(i).debug = debug_output{i};
    
    Results(i).AR1 = str2double(temp(strfind(temp,'AR1_')+4 : strfind(temp,'_AR2')-1));
    Results(i).AR2 = str2double(temp(strfind(temp,'AR2_')+4 : strfind(temp,'_tail')-1));
    Results(i).amp = str2double(temp(strfind(temp,'amp_')+4 : strfind(temp,'_lambda')-1));
    Results(i).lambda = str2double(temp(strfind(temp,'lambda_')+7 : strfind(temp,'_nlambda')-1));
    Results(i).nlambda = str2double(temp(strfind(temp,'nlambda_')+8 : strfind(temp,'_motorBC')-1));
    
end

%%
for i = 1:length(files)
    i/length(files)
    
    m = matfile([folder,files{i}],'Writable',true);
    clear results
    results.fit = Fits(i);
    results.Avg_Omega = Avg_Omega(i);
    results.Motor_Torque = Motor_Torque(i);
    
    m.results = results;
    
end
