




if actively_flipping_tail
    temp = 0; % keep final step e.g. 128 in since it isn't the same as first step
else
    temp = 1; % remove final step since it's the same as the first
end

step_list = 0:(n_phase_pts - temp);

total_period = 1 / 46.0 ;
time_list = linspace(0, total_period, n_phase_pts + 1)  ;
time_list = time_list(1:end - temp);  % leave off last value which is same as first

phase_list = linspace(0, 2*pi, n_phase_pts + 1);   % *always* 2 pi rad in a cycle
phase_list = phase_list(1:end - temp);  % leave off last value which is same as first


names = {'Body', 'Transverse', 'Tail', 'Coplanar_Hairs' , 'Normal_Top_Hairs' , 'Normal_Bottom_Hairs'};



do_figure = false;
do_vid = false;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Inds Times Files Steps
name = 'Metadata';
files = dir([folder_processed,name,'*']);
files = {files.name};

clear times steps phases
for i = 1:length(files)
    ind0 = strfind(files{i},'step_') + 5;
    ind01 = strfind(files{i},'_time') - 1;
    ind1 = strfind(files{i},'time_') + 5;
    ind2 = strfind(files{i},'_phase') - 1;
    ind3 = strfind(files{i},'.mat') - 1;
    times(i) = str2double(files{i}(ind1:ind2));
    steps(i) = str2double(files{i}(ind0:ind01));
    phases(i) = str2double(files{i}(ind2+8:ind3));
end

[~,inds] = sort(times);

Files = files(inds);
Times = times(inds);
Steps = steps(inds);
Phases = phases(inds);


% although it seems like this wouldn't work for actively_flipping_tail,
% because we flip the tail here and not in Salome and otherwise things are
% still periodic, this should work even for this case
clear Steps_closest
Steps_periodic = [Steps (Steps + n_phase_pts) ]; % create ficticious periodic version of steps so we can wrap around to beginning steps if those are closest
steps_missing = setdiff(step_list, Steps);
for i = 1:length(steps_missing)
    [inds] = find( min(abs( steps_missing(i) - Steps_periodic)) == abs( steps_missing(i) - Steps_periodic) );
    if length(inds) > 1
        inds = inds(randperm(length(inds)));  % randomize which existing step is chosen if missing step is equally close to two existing steps
    end
    Steps_closest(i) = Steps_periodic(inds(1));
end
Steps_closest(Steps_closest > (n_phase_pts - 1) ) = Steps_closest(Steps_closest > (n_phase_pts - 1) )  -  n_phase_pts;  % convert back from ficticious periodic steps to actual steps


%%  get BCs for tail and then transverse for one phase angle at a time
c = 0;
for f =  steps_missing(1:end)
    c = c + 1;
    
    
    disp(['Missing Step ',num2str(c),' out of ',num2str(length(steps_missing))]);
    
    
    step_closest = Steps_closest(c);
    metafile = Files{Steps == step_closest};
    %     load([folder, metafile ]);
    
    time = time_list(step_list == f);  %time going with this missing step
    phase = phase_list(step_list == f);
    
    
     parameters_file = [folder_original , 'Parameters_step_',num2str(step_closest), '_time_',num2str(Times(Steps == step_closest),'%.15f'),'_phase_',num2str(Phases(Steps == step_closest),'%.15f'),'.txt'];
    
%     parameters_file = 'C:/Users/rudi/Desktop/RD/dino_mesh_automation/phases_hairs_4_1/original/Parameters_step_0_time_0.000000000000000_phase_0.000000000000000.txt';
    if ~exist(parameters_file,'file')
        disp('parameters_file for this step doesn''t exist; using any parameters file found');
        temp = dir([folder_original,'Parameters_step_*.txt']);  temp = {temp.name};
        parameters_file = [folder_original, temp{1}];
    end
    dino_geom_parameters_from_file;
    
    
    clear Mesh Metadata
    
    load([folder_processed, metafile ]);
    
    for n = 1:length(names)
        meshname = [names{n}, metafile(9:end-3), 'dat'];
        
        [Mesh(n)] = load_mesh([folder_processed,meshname],[],[],'mesh');  %save rand vertex permutation for use later, so all future calculations use same exact mesh!
    end
    
    for i = 1:length(names)
        Mesh(i).name = names{i};
        
        Mesh(i).orientation = [1 0 0; 0 1 0; 0 0 1;]';
        Mesh(i).refpoints = [0 0 0]';
    end
    
    
    Mesh_orig = Mesh;
    
    
    
    [Mesh] = shift_mesh_indices(names,Mesh);  
    
    
    
    [Mesh] = global_inds(Mesh);
    [Mesh] = renumber_Mesh(Mesh);
    
    if store_constants
        clear temp_input
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
    
    
    
    
    
    %% tail
    ni = find(strcmp(names,'Tail'));
    % parameters_tail = Metadata.tail.parameters;
    
    [~,inds] = ismember(Mesh(ni).indices.orig.vert,  Metadata.Tail.indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs in Metadata to match current order in Mesh
    parameters_tail = Metadata.Tail.parameters(inds,:);
    
    
    
    verts = NaN(size(parameters_tail,1),3);  vels = verts;
    parfor p = 1:size(parameters_tail,1)
%         p/size(parameters_tail,1)
        if static_tail
            [pt, vel, der, V_u, radius] = tail_parameterized(parameters_tail(p,1), time_tail, geom.tail);
        else
            [pt, vel, der, V_u, radius] = tail_parameterized(parameters_tail(p,1), time, geom.tail);
        end
        
        pts(p,:) = pt;
        
        if ~tail_as_sheet
            radiuses(p) = radius;
            
            V_u = V_u / sqrt(sum(V_u.^2));
            der = der / sqrt(sum(der.^2));
            
            ders(p,:) = der;
            V_us(p,:) = V_u;
            
            V_surface = rotate_arbitrary_axis(V_u, [0 0 0]', der, parameters_tail(p,2));
            V_surface = V_surface / sqrt(sum(V_surface.^2));
            
            verts(p,:) = pt(:) + radius *  V_surface(:) ;
        else
            radius = parameters_tail(p,2);
            radiuses(p) = radius;
            verts(p,:) = pt(:) + [0 radius 0]'  ;  %vert on sheet must have same x and z coord as closest centerline pt
            
        end
        
        if flipped_tail
            verts(p,:) = rotate_arbitrary_axis( verts(p,:), geom.tail.translation, [0 1 0]', pi/4);
        end
        
        if actively_flipping_tail
             verts(p,:) = rotate_arbitrary_axis( verts(p,:), geom.tail.translation, [0 1 0]', slmeval(phase, tail_spline));
            
        end
        
        
        if lateral_tail
            verts(p,:) = rotate_arbitrary_axis(verts(p,:), geom.tail.translation, [0 0 1]', -25*pi/180);
        end
        if centered_tail
            verts(p,:)  =   verts(p,:)  + geom.tail.shift';
        end
        
         verts(p,:) = verts(p,:) + [-tail_shift 0 0];
         
        if static_tail
            vels(p,:) = [0 0 0];
        else
            if ~tail_as_sheet
                [n_t] = surface_vector_vel(parameters_tail(p,1), parameters_tail(p,2), time, geom.tail);  %time derivative of direction vector from centerline point to surface point
                
                vels(p,:) = vel + radius * n_t;
            else
                vels(p,:) = vel;  %velocity is constant across sheet
            end
            
        end
        
        if flipped_tail
            vels(p,:) = rotate_arbitrary_axis(vels(p,:),[0 0 0]', [0 1 0]', pi/4);  %rotate around origin since this is a velocity vector
        end
        
        if actively_flipping_tail
             vels(p,:) = rotate_arbitrary_axis(vels(p,:),[0 0 0]', [0 1 0]',  slmeval(phase, tail_spline));  %rotate around origin since this is a velocity vector
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
    
    
    
    [is_intersected] = submesh_intersections(Mesh, geom.tail.radius / 8 ,[], true, input.performance.nthreads);
    if is_intersected
        error('Self-intersections between submeshes detected');
    end
    
    
    %% transverse
    ni = find(strcmp(names,'Transverse'));
    %     parameters_transverse = Metadata.transverse.parameters;
    [~,inds] = ismember(Mesh(ni).indices.orig.vert,  Metadata.Transverse.indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs in Metadata to match current order in Mesh
    parameters_transverse = Metadata.Transverse.parameters(inds,:);
    
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
    clear hair_vels hair_parameters
    hair_cases = {'Coplanar_Hairs','Normal_Top_Hairs','Normal_Bottom_Hairs'};
    for hc = 1:length(hair_cases)
        hair_case = hair_cases{hc};
        
        
        ni = find(strcmp(names,hair_case));
        %     parameters_wingtip = Metadata.wingtip.parameters;
        [~,inds] = ismember(Mesh(ni).indices.orig.vert,  Metadata.(hair_case).indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs in Metadata to match current order in Mesh
        parameters_wingtip = Metadata.(hair_case).parameters(inds,:);
        
        verts = NaN(size(parameters_wingtip,1),3);  vels = verts;
        parfor p = 1:size(parameters_wingtip,1)
            
            [verts(p,:), vels(p,:) ] = transverse_hairs_parameterized(parameters_wingtip(p,1), parameters_wingtip(p,2), time, geom.transverse, hair_case )
            
        end
        
        Mesh(ni).verts = verts;
        wingtip_vel = vels;
        
        % make sure shared verts that were altered match exactly across all submeshes
        for nii = setdiff(1:length(Mesh),ni)
            [is_shared, inds] = ismember( Mesh(nii).indices.orig.vert , Mesh(ni).indices.orig.vert );
            
            Mesh(nii).verts( is_shared  , :) = Mesh(ni).verts( inds(is_shared), :);
        end
        
        hair_vels.(hair_case) = wingtip_vel;
        hair_parameters.(hair_case) = parameters_wingtip;
        
    end  % hair cases loop
    
    
    
    %%
    basename = ['_step_',num2str(f),'_time_', sprintf('%0.15f',time),'_phase_',sprintf('%0.15f',phase)];
    
    
    for i = 1:length(names)
        
        % delete old failed mesh file, if exists
        temp  = dir([folder_interpolated,names{i},'_step_',num2str(f),'_*.dat']);
        if ~isempty(temp)
            temp = temp.name;
            delete([folder_interpolated,temp]);
        end
        
        filename = [folder_interpolated, names{i},basename,'.dat'];
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
    
    for i = 1:length(names)  %ensures order within Metadata matches order in which submeshes were loaded (which is order of names_L)
        name = names{i};
        switch name
            case 'Body'
                Metadata.Body.BCs = zeros(Mesh(i).n_vert,3);
                
            case 'Transverse'
                Metadata.Transverse.BCs = transverse_vel;
                Metadata.Transverse.parameters = parameters_transverse;
            case 'Tail'
                Metadata.Tail.BCs = tail_vel;
                Metadata.Tail.parameters = parameters_tail;
            case 'Coplanar_Hairs'
                Metadata.Coplanar_Hairs.BCs = hair_vels.Coplanar_Hairs;
                Metadata.Coplanar_Hairs.parameters = hair_parameters.Coplanar_Hairs;
            case 'Normal_Top_Hairs'
                Metadata.Normal_Top_Hairs.BCs = hair_vels.Normal_Top_Hairs;
                Metadata.Normal_Top_Hairs.parameters = hair_parameters.Normal_Top_Hairs;
            case 'Normal_Bottom_Hairs'
                Metadata.Normal_Bottom_Hairs.BCs = hair_vels.Normal_Bottom_Hairs;
                Metadata.Normal_Bottom_Hairs.parameters = hair_parameters.Normal_Bottom_Hairs;
                
        end
        Metadata.(name).indices.orig.vert = Mesh(i).indices.orig.vert;
        %         Metadata.(name).rand_inds = Metadata_temp(i).mesh.rand_inds;
        %  Metadata.(name).rand_inds = [1:Mesh(i).n_vert];
    end
    
    filename = [folder_interpolated, 'Metadata',basename,'.mat'];
    
    save(filename,'Metadata');
    
    
    
end

