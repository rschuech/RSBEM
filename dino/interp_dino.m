function [Mesh, Vels, Parameters, Metadata] = interp_dino(time, folder, metafile, subnames, geom, nthreads, check_intersections , varargin)

% names_L = lower(subnames);
names_L = subnames;

store_constants = true;

load([folder, metafile ]);  % loads Metadata

if ~isempty(varargin)
    rand_inds = varargin{1};
    names = {rand_inds.name};
end

for n = 1:length(subnames)
    meshname = [subnames{n}, metafile(9:end-3), 'dat'];
    if isempty(varargin)
        [Mesh(n)] = load_mesh([folder,meshname],[],[],'mesh');  %save rand vertex permutation for use later, so all future calculations use same exact mesh!
    else
        rand_inds_ind = find(strcmp(names_L{n},names));
        [Mesh(n)] = load_mesh([folder,meshname],[],[rand_inds(rand_inds_ind).mesh.rand_inds],'mesh');
    end
end

for i = 1:length(subnames)
    Mesh(i).name = names_L{i};
    
    Mesh(i).orientation = [1 0 0; 0 1 0; 0 0 1;]';
    Mesh(i).refpoints = [0 0 0]';
end

Mesh_orig = Mesh;

ind = find(strcmp('Tail',names_L));  %Mesh ind going with tail, if we loaded tail
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
    
    temp_input.performance.nthreads = nthreads;
    
    temp_input.accuracy.integration_tol.area.abstol = 1E-9 ;
    temp_input.accuracy.integration_tol.area.reltol = 0.1;
    temp_input.accuracy.integration_tol.area.maxevals = Inf;
    
    temp_input.accuracy.integration_tol.centroid.abstol = 1E-9;
    temp_input.accuracy.integration_tol.centroid.reltol = 0.1;
    temp_input.accuracy.integration_tol.centroid.maxevals = Inf;
    
    temp_input.accuracy.integration_tol.volume.abstol = 1E-9;
    temp_input.accuracy.integration_tol.volume.reltol = 0.1;
    temp_input.accuracy.integration_tol.volume.maxevals = Inf;
    
    [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
end

if ismember('tail',names_L) && check_intersections
    [is_intersected] = submesh_intersections(Mesh, geom.tail.radius / 8 ,[], true, nthreads);
    if is_intersected
        stopafra
    end
end


%% tail
ni = find(strcmp(subnames,'Tail'));
% parameters_tail = Metadata.tail.parameters;

[~,inds] = ismember(Mesh(ni).indices.orig.vert,  Metadata.Tail.indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs in Metadata to match current order in Mesh
parameters_tail = Metadata.Tail.parameters(inds,:);



verts = NaN(size(parameters_tail,1),3);  vels = verts;
parfor p = 1:size(parameters_tail,1)
    [pt, vel, der, V_u, radius] = tail_parameterized(parameters_tail(p,1), time, geom.tail);
    
    pts(p,:) = pt;
    
    V_u = V_u / sqrt(sum(V_u.^2));
    der = der / sqrt(sum(der.^2));
    
    ders(p,:) = der;
    V_us(p,:) = V_u;
    
    V_surface = rotate_arbitrary_axis(V_u, [0 0 0]', der, parameters_tail(p,2));
    V_surface = V_surface / sqrt(sum(V_surface.^2));
    
    verts(p,:) = pt(:) + radius *  V_surface(:) ;
    
    %         rotated = rotate_arbitrary_axis(input, point, vec, angle)
    
    [n_t] = surface_vector_vel(parameters_tail(p,1), parameters_tail(p,2), time, geom.tail);  %time derivative of direction vector from centerline point to surface point
    
    vels(p,:) = vel + radius * n_t;
   
    
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
ni = find(strcmp(subnames,'Transverse'));
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

%% hair sheets
clear hair_vels hair_parameters
hair_cases = {'Coplanar_Hairs','Normal_Top_Hairs','Normal_Bottom_Hairs'};
for hc = 1:length(hair_cases)
    hair_case = hair_cases{hc};
    
    
    ni = find(strcmp(subnames,hair_case));
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
    
    Vels.(hair_case) = wingtip_vel;
    Parameters.(hair_case) = parameters_wingtip;
    
end


Vels.tail = tail_vel;
Vels.transverse = transverse_vel;
%     Vels.wingtip = wingtip_vel;

Parameters.tail = parameters_tail;
Parameters.transverse = parameters_transverse;
% Parameters.wingtip = parameters_wingtip;