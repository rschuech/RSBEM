load('C:\Users\rudi\Desktop\RD\bad_inputs.mat')


folder = 'C:\Users\rudi\Desktop\RD\minidumps_archived\';
folder = 'C:\Users\rudi\Desktop\RD\swept_dumps_archived\';

for i = 1:length(bad_inputs)
    
    files = dir([folder,'*',bad_inputs(i).paths.namebase.body,'_*',bad_inputs(i).paths.namebase.tail,'_*']);
    
    files = {files.name};
    
    for j = 1:length(files)
        [folder,files{j}]
       %  pause
        delete([folder,files{j}]);
    end
    
    
end



%% check meshes

folder = 'E:\Hull\swept_meshes\';
folder = 'C:\Users\rudi\Desktop\RD\swept_meshes\';

files = dir([folder,'*metadata*']);

files = {files.name};

for i = 1:length(files)
    temp = load([folder,files{i}]);
    if isfield(temp.Metadata.mesh,'meshing_succeeded') &&   ~temp.Metadata.mesh.meshing_succeeded
    files{i}
    delete([folder,files{i}]);
    end
end
