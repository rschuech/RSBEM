
folder = 'C:\Users\rudi\Desktop\RD\swept_dumps\';

% files = dir([folder,'*time_dump.mat']);
% files = {files.name};
%
% F = length(files);


for f = 1:length(Results)  %timestepping dumps
    forced_file = [Results(f).name,'_forced_dump.mat'];
    f/length(Results)
    
    
    if exist([folder,forced_file])
        temp = load([folder,forced_file],'fcoeffs');
        
        Results(f).fcoeffs.translation = temp.fcoeffs.translation;
        Results(f).fcoeffs.rotation = temp.fcoeffs.rotation;
    else
        Results(f).fcoeffs.translation = NaN(1,3);
        Results(f).fcoeffs.rotation = NaN(1,3);
    end
    
end



