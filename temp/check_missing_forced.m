% folder = 'C:\Users\rudi\Desktop\RD\swept_dumps_archived\';
folder = 'C:\Users\rudi\Desktop\RD\swept_dumps\';
% folder2 = 'E:\Hull\CFD03\swept_dumps\';



files = dir([folder,'*torque_dump.mat']);
files = {files.name};

% clear Fits Avg_Omega Motor_Torque
%   Fits = struct('T',cell(1,length(files)),'line',cell(1,length(files)),'converged',cell(1,length(files)));
F = length(files);
% Results.file = cell(F,1);

% files = files(randperm(length(files)));

missing_forced = false(size(files));

% ppm = ParforProgressStarter2('shat', length(files));

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
    
    forced_file = [dump_file(1:end-23), 'forced_dump.mat'];
    if ~exist([folder, forced_file]) %no partial file, load timestepping file if it exists
       missing_forced(f) = true;
    end
    
end