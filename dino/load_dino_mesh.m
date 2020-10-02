function [Mesh, Metadata, matrix_props , Metadata_rand_inds] = load_dino_mesh(phase, Mesh_files, input, varargin)

places = 13;  %how many decimal places to round phase angles to when comparing - necessary due to limited precision of filenames

ind = find(roundn(phase,-places) == roundn(Mesh_files.phase,-places));
if isempty(ind)
    disp('Can''t find mesh with same phase!');
    Mesh = [];  matrix_props = [];  Metadata_rand_inds= [];  original_indices = [];
    return
elseif length(ind) > 1
    disp('Multiple meshes found with same phase!');
    Mesh = [];  matrix_props = [];  Metadata_rand_inds= []; original_indices = [];
    return
end

load([input.paths.datfolder,Mesh_files.Metadata{ind}]);  %this sneakily loads Metadata

if ~isfield(Metadata.geom,'phase_speed')
    Metadata.geom.phase_speed =   2 * pi * 46;
    disp('hardcoding phase speed in load_dino_mesh.m')
end

switch input.bugtype
    case 'dino'
        
        fields = setdiff( fieldnames(Mesh_files) , {'phase','time','Metadata'});
        for i = 1:length(fields)
            field = fields{i};
            if isempty(varargin)
%                 [Mesh(i), Metadata_rand_inds(i)] = load_mesh([input.paths.datfolder,Mesh_files.(field){ind}],[],[],'mesh');
            [Mesh(i), Metadata_rand_inds(i), index_mapping.local2global{i}] = load_mesh([input.paths.datfolder,Mesh_files.(field){ind}],[],[],'mesh');
            else
                temp = varargin{1};
                rand_inds_names = {temp.name};
                rand_inds = temp(strcmp(rand_inds_names, field)).mesh.rand_inds;
%                 [Mesh(i), Metadata_rand_inds(i)] = load_mesh([input.paths.datfolder,Mesh_files.(field){ind}],[],rand_inds,'mesh');
                 [Mesh(i), Metadata_rand_inds(i), index_mapping.local2global{i}] = load_mesh([input.paths.datfolder,Mesh_files.(field){ind}],[],rand_inds,'mesh');
            end
        end
      
        for i = 1:length(fields)
            Mesh(i).name = fields{i};
            Mesh(i).parent_topology = input.parent_topology.( Mesh(i).name );
            Metadata_rand_inds(i).name = fields{i};
        end
        

     % eliminate any possible overlap in global indices across coincident submesh groups, and also make sure global indices are consecutive
    [Mesh, index_mapping.local2global] = shift_global_indices(Mesh, index_mapping.local2global, input.coincident_submeshes);
    
        
        %since Metadata is always created by loading all 3 submeshes, need to match
        %that convention here (i.e. tail submesh indices are always shifted up
        %using both body and transverse indices)
        
        temp = {'Body','Transverse','Tail','Coplanar_Hairs','Normal_Top_Hairs','Normal_Bottom_Hairs'};
        removes = temp( ~logical(input.potatohead));
        [~,inds] = ismember(removes, fields);  % which submesh inds (that were loaded at all in first place) do we have to remove?
        Mesh(inds(inds ~= 0)) = [];
      
        
    case 'sheet'
        Mesh(1) = load_mesh([input.paths.datfolder,Mesh_files.sheet{ind}],[],[],'mesh');
        Mesh(1).name = 'sheet';
        
end


switch input.bugtype
    case 'dino'
        
        for si = 1:length(Mesh)
            Mesh(si).orientation = [1 0 0; 0 1 0; 0 0 1;]';
            % set all refpoints the same for ease of comparison
%             if isfield(Metadata.geom,'body')
                Mesh(si).refpoints = [Metadata.geom.body.center];
%             else
%                 Mesh(si).refpoints = [7 0 0]';  %more or less center of transverse
%             end
        end
        
    case 'sheet'
        Mesh.orientation = [1 0 0; 0 1 0; 0 0 1;]';
        % set all refpoints the same for ease of comparison
       if isfield(Metadata.geom.sheet,'center')
           Mesh.refpoints = [Metadata.geom.sheet.center];
       else
           Mesh.refpoints = [30 0 0]';
       end
end

[Mesh] = global_inds(Mesh);
[Mesh] = renumber_Mesh(Mesh);


temp_input.performance.nthreads = input.performance.nthreads;
temp_input.accuracy.integration_tol.area.abstol = input.accuracy.integration_tol.area.abstol ;
temp_input.accuracy.integration_tol.area.reltol = input.accuracy.integration_tol.area.reltol ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
temp_input.accuracy.integration_tol.area.maxevals = input.accuracy.integration_tol.area.maxevals;
temp_input.accuracy.integration_tol.centroid.abstol = input.accuracy.integration_tol.centroid.abstol ;
temp_input.accuracy.integration_tol.centroid.reltol = input.accuracy.integration_tol.centroid.reltol; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
temp_input.accuracy.integration_tol.centroid.maxevals = input.accuracy.integration_tol.centroid.maxevals;
temp_input.accuracy.integration_tol.volume.abstol = input.accuracy.integration_tol.volume.abstol ;
temp_input.accuracy.integration_tol.volume.reltol = input.accuracy.integration_tol.volume.reltol ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
temp_input.accuracy.integration_tol.volume.maxevals = input.accuracy.integration_tol.volume.maxevals;

[Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);

[matrix_props] = gen_matrix_props(input,Mesh);