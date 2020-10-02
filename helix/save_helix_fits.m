folder = 'E:\Hull\CFD03\swept_dumps\';
folder2 = 'E:\Hull\CFD03\swept_dumps\';

files = dir([folder,'*timestepping.mat']);
files = {files.name};
clear Fits

parfor f = 1:length(files)
    file = files{f};
    f/length(files)
    
    newfile = [file(1:end-16) 'helix_fit.mat'];
    %
    %     if exist([folder2,newfile],'file')
    %         continue
    %     end
    
    temp = load([folder,file]);
    
    timestepping_solution = temp.timestepping_solution;
    
    [Fits{f}] = helix_fitting(timestepping_solution);
    
   % save([folder2,newfile],'fits');
    
end