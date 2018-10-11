function [Mesh_files] = get_sheet_mesh_files(folder)

%folder = input.paths.datfolder;

clear Inds Files Times Phases
names = {'Sheet','Metadata'};
for n = 1:length(names)
    name = names{n};
    
    files = dir([folder,name,'*']);
    files = {files.name};
    %files = files(3:end);
    
    clear times phases
    for i = 1:length(files)
        ind1 = strfind(files{i},'time_') + 5;
        ind2 = strfind(files{i},'_phase') - 1;
        ind3 = strfind(files{i},'.dat') - 1;
        times(i) = str2double(files{i}(ind1:ind2));
        phases(i) = str2double(files{i}(ind2+8:ind3));
    end
    
    [~,inds] = sort(times);
    Inds{n} = inds;
    Files{n} = files(inds);
    Times{n} = times(inds);
    Phases{n} = phases(inds);
end


    Time = Times{1};
    Phase = Phases{1};


Mesh_files.time = Time;
Mesh_files.phase = Phase;
Mesh_files.sheet = Files{1};
Mesh_files.metadata = Files{2};
%%
