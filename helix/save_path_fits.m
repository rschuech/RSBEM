folder = 'E:\Hull\CFD03\swept_dumps\';
folder2 = 'E:\Hull\CFD03\swept_dumps\';

files = dir([folder,'*timestepping.mat']);
files = {files.name};

clear Fits Avg_Omega Motor_Torque
%   Fits = struct('T',cell(1,length(files)),'line',cell(1,length(files)),'converged',cell(1,length(files)));

ppm = ParforProgressStarter2('shat', length(files));

parfor f = 1:length(files)  %timestepping dumps
    timestepping_file = files{f};
    f/length(files)
    
    interp_file = [ timestepping_file(1:end-16) 'dump.mat'];
    
    if ~exist([folder,interp_file])
        continue
    end
    
    newfile = [timestepping_file(1:end-16) 'helix_fit.mat'];
    %
    %     if exist([folder2,newfile],'file')
    %         continue
    %     end
    
    temp = load([folder,timestepping_file]);
    
    timestepping_solution = temp.timestepping_solution;
    
    [temp] = path_fitting_old(timestepping_solution);
   fields = fieldnames(temp);
   for ff = 1:length(fields)
    Fits(f).(fields{ff}) = temp.(fields{ff});
   end
 
    
    %% get avg omega from interpolant in interpolant dump file
    temp = load([folder,interp_file],'interpolant','input');
    interpolant = temp.interpolant(end).omega;  %last and most accurate interpolant for tail rotation rate omega
    
    [avg_omega] = compute_avg_omega(interpolant);
    
    Avg_Omega(f) = avg_omega;
    
    Motor_Torque(f) = temp.input.tail.motor_torque;
    % save([folder2,newfile],'fits');
    ppm.increment(f);
    
    
end

delete(ppm);