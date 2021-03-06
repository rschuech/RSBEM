
folder = 'E:\Hull\dinoflagellate\meshes_transverse_hairs_parallel2_3_all\';
folder = 'E:\Hull\dinoflagellate\meshes_parallel2_3\';

folder = 'E:\Hull\dinoflagellate\meshes_parallel2_3_modded\';
outfolder = 'E:\Hull\dinoflagellate\meshes_parallel2_3_modded_interped\';

folder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_3_modded\';
outfolder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_3_all3\';

folder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_3_filtered\';
outfolder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_3_filtered_BCs\';

folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_modded\';
outfolder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_modded\';

folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_interpBCs\';
outfolder = folder;

folder = 'C:\Users\rudi\Desktop\RD\meshes_deep_groove_BCs\';
outfolder = folder;


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


static_tail = true;
 time_tail = 0;  %must match tail mesh file we are using for all time
flipped_tail = false;

static_tail = false;
flipped_tail = false;
lateral_tail = false;


input.performance.nthreads = 40;
store_constants = true;


dino_geom_parameters_huge_groove;

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
for f = 0% steps_missing(1:end)
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
    
    
    
    
    
    
    
    
    
clear Mesh Metadata

load([folder, metafile ]);

for n = 1:length(names)
    meshname = [names{n}, metafile(9:end-3), 'dat'];
    
    [Mesh(n)] = load_mesh([folder,meshname],[],[],'mesh');  %save rand vertex permutation for use later, so all future calculations use same exact mesh!
end

for i = 1:length(names)
    Mesh(i).name = names_L{i};
    
    Mesh(i).orientation = [1 0 0; 0 1 0; 0 0 1;]';
    Mesh(i).refpoints = [0 0 0]';
end


Mesh_orig = Mesh;

    
    
    ind = find(strcmp('tail',names_L));  %Mesh ind going with tail, if we loaded tail
    if ~isempty(ind)
        other_inds = setdiff(1:length(names_L), ind);  % inds for all other submeshes
        %tail indices might overlap body/transverse/wingtip since they're not part of same
        %mesh
        other_orig_elem = [];  other_orig_vert = [];
        for other_ind = other_inds
            other_orig_elem = [other_orig_elem; Mesh(other_ind).indices.orig.elem];
            other_orig_vert = [other_orig_vert; Mesh(other_ind).indices.orig.vert];
        end
        
        
        Mesh(ind).indices.orig.elem = Mesh(ind).indices.orig.elem + max(other_orig_elem);
        Mesh(ind).indices.orig.vert = Mesh(ind).indices.orig.vert + max(other_orig_vert);
        Mesh(ind).elems             = Mesh(ind).elems             + max(other_orig_vert);  %still using orig vert labels here, until renumber_Mesh below
    end
    

    
    
    [Mesh] = global_inds(Mesh);
    [Mesh] = renumber_Mesh(Mesh);
    
    if store_constants
        
        temp_input.performance.nthreads = input.performance.nthreads;
        
        temp_input.accuracy.integration_tol.area.abstol = 1E-9 ;
        temp_input.accuracy.integration_tol.area.reltol = 0.1;
        temp_input.accuracy.integration_tol.area.maxevals = Inf;  
        
        temp_input.accuracy.integration_tol.centroid.abstol = 1E-9;
        temp_input.accuracy.integration_tol.centroid.reltol = 0.1 ;
        temp_input.accuracy.integration_tol.centroid.maxevals = Inf;  
        
        temp_input.accuracy.integration_tol.volume.abstol = 1E-9;
        temp_input.accuracy.integration_tol.volume.reltol = 0.1 ;
        temp_input.accuracy.integration_tol.volume.maxevals = Inf; 
        
        [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
    end
    
    if ismember('tail',names_L)
        [is_intersected] = submesh_intersections(Mesh, geom.tail.radius / 8 ,[], true, input.performance.nthreads);
        if is_intersected
            stopafra
        end
    end
    
    
    
    
    
    
    
    %% tail
    ni = find(strcmp(names,'Tail'));
    % parameters_tail = Metadata.tail.parameters;
    
    [~,inds] = ismember(Mesh(ni).indices.orig.vert,  Metadata.tail.indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs in Metadata to match current order in Mesh
    parameters_tail = Metadata.tail.parameters(inds,:);
    
    
    
    verts = NaN(size(parameters_tail,1),3);  vels = verts;
    parfor p = 1:size(parameters_tail,1)
        if static_tail
        [pt, vel, der, V_u, radius] = tail_parameterized(parameters_tail(p,1), time_tail, geom.tail);
        else
            [pt, vel, der, V_u, radius] = tail_parameterized(parameters_tail(p,1), time, geom.tail);
        end
        radiuses(p) = radius;
        pts(p,:) = pt;
        
        V_u = V_u / sqrt(sum(V_u.^2));
        der = der / sqrt(sum(der.^2));
        
        ders(p,:) = der;
        V_us(p,:) = V_u;
        
        V_surface = rotate_arbitrary_axis(V_u, [0 0 0]', der, parameters_tail(p,2));
        V_surface = V_surface / sqrt(sum(V_surface.^2));
        
        verts(p,:) = pt(:) + radius *  V_surface(:) ;
        
          if flipped_tail
                verts(p,:) = rotate_arbitrary_axis( verts(p,:), geom.tail.translation, [0 1 0]', pi/2);
          end
                 if lateral_tail
                verts(p,:) = rotate_arbitrary_axis(verts(p,:), geom.tail.translation, [0 0 1]', -25*pi/180);
           end
       
           
        if static_tail
            vels(p,:) = [0 0 0];
        else
            
            %         rotated = rotate_arbitrary_axis(input, point, vec, angle)
            
            [n_t] = surface_vector_vel(parameters_tail(p,1), parameters_tail(p,2), time, geom.tail);  %time derivative of direction vector from centerline point to surface point
            
            vels(p,:) = vel + radius * n_t;
            
            
        end
        
        if flipped_tail
            vels(p,:) = rotate_arbitrary_axis(vels(p,:),[0 0 0]', [0 1 0]', pi/2);  %rotate around origin since this is a velocity vector
        end
                 if lateral_tail  % using 25 degrees since measured was 24 degrees
                vels(p,:) = rotate_arbitrary_axis(vels(p,:),[0 0 0]', [0 0 1]', -25*pi/180);  %rotate around origin since this is a velocity vector
            end
        
        
    end
    % sometimes there is apparently a very small imaginary component that
    % gets introduced by above operations and this makes mex angry
    Mesh(ni).verts = real(verts);
    tail_vel = real(vels);
    
    % make sure shared verts that were altered match exactly across all submeshes
    for nii = setdiff(1:length(Mesh),ni)
        [is_shared, inds] = ismember( Mesh(nii).indices.orig.vert , Mesh(ni).indices.orig.vert );
        
        Mesh(nii).verts( is_shared  , :) = Mesh(ni).verts( inds(is_shared), :);
    end
    
    %% transverse
    ni = find(strcmp(names,'Transverse'));
    %     parameters_transverse = Metadata.transverse.parameters;
    [~,inds] = ismember(Mesh(ni).indices.orig.vert,  Metadata.transverse.indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs in Metadata to match current order in Mesh
    parameters_transverse = Metadata.transverse.parameters(inds,:);
    
    verts = NaN(size(parameters_transverse,1),3);  vels = verts;
    parfor p = 1:size(parameters_transverse,1)
        
        [verts(p,:), vels(p,:) ] = transverse_parameterized(parameters_transverse(p,1),parameters_transverse(p,2), time, geom.transverse)
        
        
        
    end
    
    Mesh(ni).verts = verts;
    transverse_vel = vels;
    
    % make sure shared verts that were altered match exactly across all submeshes
    for nii = setdiff(1:length(Mesh),ni)
        [is_shared, inds] = ismember( Mesh(nii).indices.orig.vert , Mesh(ni).indices.orig.vert );
        
        Mesh(nii).verts( is_shared  , :) = Mesh(ni).verts( inds(is_shared), :);
    end
    
    %% wingtip
    ni = find(strcmp(names,'Wingtip'));
%     parameters_wingtip = Metadata.wingtip.parameters;
             [~,inds] = ismember(Mesh(ni).indices.orig.vert,  Metadata.wingtip.indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs in Metadata to match current order in Mesh
                    parameters_wingtip = Metadata.wingtip.parameters(inds,:);
                    
    verts = NaN(size(parameters_wingtip,1),3);  vels = verts;
    parfor p = 1:size(parameters_wingtip,1)
        
        [verts(p,:), vels(p,:) ] = transverse_hairs_parameterized(parameters_wingtip(p,1), parameters_wingtip(p,2), time, geom.transverse )
        
        
        
    end
    
    Mesh(ni).verts = verts;
    wingtip_vel = vels;
    
    % make sure shared verts that were altered match exactly across all submeshes
    for nii = setdiff(1:length(Mesh),ni)
        [is_shared, inds] = ismember( Mesh(nii).indices.orig.vert , Mesh(ni).indices.orig.vert );
        
        Mesh(nii).verts( is_shared  , :) = Mesh(ni).verts( inds(is_shared), :);
    end
    
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
                Metadata.transverse.BCs = transverse_vel;
                Metadata.transverse.parameters = parameters_transverse;
            case 'tail'
                Metadata.tail.BCs = tail_vel;
                Metadata.tail.parameters = parameters_tail;
            case 'wingtip'
                Metadata.wingtip.BCs = wingtip_vel;
                Metadata.wingtip.parameters = parameters_wingtip;
        end
        Metadata.(name).indices.orig.vert = Mesh(i).indices.orig.vert;
        %         Metadata.(name).rand_inds = Metadata_temp(i).mesh.rand_inds;
      %  Metadata.(name).rand_inds = [1:Mesh(i).n_vert];
    end
    
    filename = [outfolder, 'Metadata',basename,'.mat'];
    
    save(filename,'Metadata');
    
    
    
end