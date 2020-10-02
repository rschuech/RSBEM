
folder = 'E:\Hull\dinoflagellate\meshes_transverse_hairs_parallel2_3_all\';
folder = 'E:\Hull\dinoflagellate\meshes_parallel2_3\';

folder = 'E:\Hull\dinoflagellate\meshes_parallel2_3_modded\';
outfolder = 'E:\Hull\dinoflagellate\meshes_parallel2_3_modded_interped\';

folder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_3_modded\';
outfolder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_3_all3\';


if ~exist(outfolder,'dir')
    mkdir(outfolder);
end


n_phase_pts = 128 ;% number of equally spaced configurations and meshes to compute over a full beat cycle of both transverse and tail
step_list = 0:(n_phase_pts-1);

total_period = 1 / 46.0 ;
time_list = linspace(0, total_period, n_phase_pts + 1)  ;
time_list = time_list(1:end - 1);  % leave off last value which is same as first

phase_list = linspace(0, 2*pi, n_phase_pts + 1);   % *always* 2 pi rad in a cycle
phase_list = phase_list(1:end-1);  % leave off last value which is same as first



% names = {'Body','Transverse','Tail'};
names = {'Body', 'Transverse', 'Tail', 'Wingtip'};
names_L = {'body', 'transverse', 'tail', 'wingtip'};

% names = { 'Transverse', 'Wingtip'};
% names_L = { 'transverse', 'wingtip'};

do_figure = false;
do_vid = false;

nthreads = 40;

dino_geom_parameters;  % creates geom structure

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Inds Times Files Steps
name = 'Metadata';
files = dir([folder,name,'*']);
files = {files.name};

clear times steps
for i = 1:length(files)
    ind0 = strfind(files{i},'step_') + 5;
    ind01 = strfind(files{i},'_time') - 1;
    ind1 = strfind(files{i},'time_') + 5;
    ind2 = strfind(files{i},'_phase') - 1;
    times(i) = str2double(files{i}(ind1:ind2));
    steps(i) = str2double(files{i}(ind0:ind01));
end

[~,inds] = sort(times);

Files = files(inds);
Times = times(inds);
Steps = steps(inds);



clear Steps_closest
steps_missing = setdiff(step_list, Steps);
for i = 1:length(steps_missing)
    [inds] = find( min(abs( steps_missing(i) - Steps)) == abs( steps_missing(i) - Steps) );
    if length(inds) > 1
        inds = inds(randperm(length(inds)));  % randomize which existing step is chosen if missing step is equally close to two existing steps
    end
    Steps_closest(i) = Steps(inds(1));
end


%%  get BCs for tail and then transverse for one phase angle at a time
c = 0;
for f =   steps_missing(1:end)
    c = c + 1;
    
    
    disp(['Missing Step ',num2str(c),' out of ',num2str(length(steps_missing))]);
    
    
    step_closest = Steps_closest(c);
    metafile = Files{Steps == step_closest};
%     load([folder, metafile ]);
    
    time = time_list(step_list == f);  %time going with this missing step
    phase = phase_list(step_list == f);
    % sloppily load closest existing mesh file for slight modification
%     clear Mesh
%     for n = 1:length(names)
%         meshname = [names{n}, metafile(9:end-3), 'dat'];
%         
%         
%         temp = dlmread([folder,meshname]);
%         
%         Mesh(n).n_vert = temp(1,1);  %# vertices, may or may not be the number that's actually used - see below
%         Mesh(n).n_elem = temp(1,2); %# elements
%         
%         Mesh(n).indices.orig.vert = temp(2:Mesh(n).n_vert+1,1);
%         Mesh(n).indices.orig.elem = temp(Mesh(n).n_vert+2:end,1);
%         
%         Mesh(n).verts = temp(2:Mesh(n).n_vert+1,2:4); %vertex coords
%         
%         Mesh(n).elems = temp(Mesh(n).n_vert+2:end,3:end);  %vertex indices that form each element
%     end
 

[Mesh, Vels, Parameters] = interp_dino(time, folder, metafile, names, geom, nthreads );
    
    %%
    basename = ['_step_',num2str(f),'_time_', sprintf('%0.15f',time),'_phase_',sprintf('%0.15f',phase)];
    
    
    for i = 1:length(names)
        
        
        filename = [outfolder, names{i},basename,'.dat'];
        if exist(filename,'file')
            delete(filename);
        end
        dlmwrite(filename,[Mesh(i).n_vert Mesh(i).n_elem],'-append','delimiter',' ','precision',16)
        dlmwrite(filename,[Mesh_orig(i).indices.orig.vert    Mesh(i).verts],'-append','delimiter',' ','precision',16)
        dlmwrite(filename,[Mesh_orig(i).indices.orig.elem      repmat(206,Mesh(i).n_elem,1)     Mesh_orig(i).elems],'-append','delimiter',' ','precision',16)
        
    end
    
    %%
    clear Metadata
    Metadata.geom = geom;
    Metadata.time = time;
    
    for i = 1:length(names_L)  %ensures order within Metadata matches order in which submeshes were loaded (which is order of names_L)
        name = names_L{i};
        switch name
            case 'body'
                Metadata.body.BCs = zeros(Mesh(i).n_vert,3);
                
            case 'transverse'
                Metadata.transverse.BCs = Vels.transverse;
                Metadata.transverse.parameters = Parameters.transverse;
            case 'tail'
                Metadata.tail.BCs = Vels.tail;
                Metadata.tail.parameters = Parameters.tail;
            case 'wingtip'
                Metadata.wingtip.BCs = Vels.wingtip;
                Metadata.wingtip.parameters = Parameters.wingtip;
        end
        Metadata.(name).indices.orig.vert = Mesh(i).indices.orig.vert;
        %         Metadata.(name).rand_inds = Metadata_temp(i).mesh.rand_inds;
      %  Metadata.(name).rand_inds = [1:Mesh(i).n_vert];
    end
    
    filename = [outfolder, 'Metadata',basename,'.mat'];
    
    save(filename,'Metadata');
    
    
    
end