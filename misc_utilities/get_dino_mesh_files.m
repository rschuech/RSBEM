function [Mesh_files] = get_dino_mesh_files(folder)

%folder = input.paths.datfolder;


%first figure out which submeshes were loaded when computing BCs, since
%we'll need to load all of them regardless of potatohead settings
temp = dir([folder,'Metadata*']);
temp = {temp.name};
temp = temp{1};  %arbitrarily pick first Metadata file, they should all have same submeshes....
temp = load([folder,temp]);
temp = temp.Metadata;
submeshes = setdiff(fieldnames(temp),{'geom','time'});  %all the submeshes that were loaded in compute BCs

clear Inds Files Times Phases
% names = {'Body','Transverse','Tail','Metadata'};
names = submeshes;  names{end+1} = 'Metadata';

for n = 1:length(names)
    name = names{n};
    name(1) = upper(name(1));  % e.g. transverse to Transverse
    
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

% if ~isequal(Times{1},Times{2},Times{3})
%     disp('problemo')
%     pause
% else
    Time = Times{1};
    Phase = Phases{1};
% end

Mesh_files.time = Time;
Mesh_files.phase = Phase;

for i = 1:length(names)
    name = names{i};
    Mesh_files.((name)) = Files{i};
end

% Mesh_files.body = Files{1};
% Mesh_files.transverse = Files{2};
% Mesh_files.tail = Files{3};
% Mesh_files.metadata = Files{4};
%%
