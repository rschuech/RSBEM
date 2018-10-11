function [Mesh, Metadata, matrix_props , Metadata_rand_inds] = load_dino_mesh(phase, Mesh_files, input, varargin)

places = 13;  %how many decimal places to round phase angles to when comparing - necessary due to limited precision of filenames

ind = find(roundn(phase,-places) == roundn(Mesh_files.phase,-places));
if isempty(ind)
    disp('Can''t find mesh with same phase!');
    Mesh = [];  matrix_props = [];  Metadata_rand_inds= [];
    return
elseif length(ind) > 1
    disp('Multiple meshes found with same phase!');
    Mesh = [];  matrix_props = [];  Metadata_rand_inds= [];
    return
end

load([input.paths.datfolder,Mesh_files.Metadata{ind}]);  %this sneakily loads Metadata

if ~isfield(Metadata.geom,'phase_speed')
    Metadata.geom.phase_speed =   2 * pi * 46;
    disp('hardcoding phase speed in load_dino_mesh.m')
end


% Mesh(1) = load_mesh([input.paths.datfolder,Mesh_files.body{ind}],[],[],'mesh');
% Mesh(2) = load_mesh([input.paths.datfolder,Mesh_files.transverse{ind}],[],[],'mesh');
% Mesh(3) = load_mesh([input.paths.datfolder,Mesh_files.tail{ind}],[],[],'mesh');
%
%  Mesh(1).name = 'body';
% Mesh(2).name = 'transverse';
%  Mesh(3).name = 'tail';
switch input.bugtype
    case 'dino'
        
        fields = setdiff( fieldnames(Mesh_files) , {'phase','time','Metadata'});
        for i = 1:length(fields)
            field = fields{i};
            if isempty(varargin)
                [Mesh(i), Metadata_rand_inds(i)] = load_mesh([input.paths.datfolder,Mesh_files.(field){ind}],[],[],'mesh');
            else
                temp = varargin{1};
                rand_inds_names = {temp.name};
                rand_inds = temp(strcmp(rand_inds_names, field)).mesh.rand_inds;
                [Mesh(i), Metadata_rand_inds(i)] = load_mesh([input.paths.datfolder,Mesh_files.(field){ind}],[],rand_inds,'mesh');
            end
        end
        
        for i = 1:length(fields)
            Mesh(i).name = fields{i};
            Metadata_rand_inds(i).name = fields{i};
        end
        
        %         Mesh(1) = load_mesh([input.paths.datfolder,Mesh_files.body{ind}],[],[],'mesh');
        %         Mesh(2) = load_mesh([input.paths.datfolder,Mesh_files.transverse{ind}],[],[],'mesh');
%         Mesh(3) = load_mesh([input.paths.datfolder,Mesh_files.tail{ind}],[],[],'mesh');
%         Mesh(1).name = 'body';
%         Mesh(2).name = 'transverse';
%         Mesh(3).name = 'tail';


        %we have a tail, and it should be the last submesh
        %tail indices might overlap body/transverse since they're not part of same
        %mesh
        
%      ind = find(strcmp('tail',fields));  %Mesh ind going with tail, if we loaded tail
%     if ~isempty(ind)
%         other_inds = setdiff(1:length(fields), ind);  % inds for all other submeshes
%         %tail indices might overlap body/transverse/wingtip since they're not part of same
%         %mesh
%         other_orig_elem = [];  other_orig_vert = [];
%         for other_ind = other_inds
%             other_orig_elem = [other_orig_elem; Mesh(other_ind).indices.orig.elem];
%             other_orig_vert = [other_orig_vert; Mesh(other_ind).indices.orig.vert];
%         end
%         
%         
%         Mesh(ind).indices.orig.elem = Mesh(ind).indices.orig.elem + max(other_orig_elem);
%         Mesh(ind).indices.orig.vert = Mesh(ind).indices.orig.vert + max(other_orig_vert);
%         Mesh(ind).elems             = Mesh(ind).elems             + max(other_orig_vert);  %still using orig vert labels here, until renumber_Mesh below
%     end
    
    [Mesh] = shift_mesh_indices(fields,Mesh);  % replaces above shat
    
%         Mesh(3).indices.orig.elem = Mesh(3).indices.orig.elem + max([Mesh(1).indices.orig.elem; Mesh(2).indices.orig.elem]);
%         Mesh(3).indices.orig.vert = Mesh(3).indices.orig.vert + max([Mesh(1).indices.orig.vert; Mesh(2).indices.orig.vert]);
%         Mesh(3).elems             = Mesh(3).elems             + max([Mesh(1).indices.orig.vert; Mesh(2).indices.orig.vert]);  %still using orig vert labels here, until renumber_Mesh below
        
        
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
%
% Mesh(1) = load_mesh([input.paths.datfolder,Mesh_files.tail{ind}],[],[],'mesh');
%  Mesh(1).name = 'tail';

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

% Mesh(3).refpoints = Metadata.geom.tail.translation; %this is tail_join_pt, which is center of starting hemisphere




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

clear matrix_props
% matrix_props.n_col = sum([Mesh.n_vert]);  %number collocation points in A matrix, could be greater than number of points at which to determine traction, in theory
%    matrix_props.n_col = Mesh(end).indices.glob.vert( find( ~isnan(Mesh(end).indices.glob.vert), 1, 'last' ) );  %highest global index must be in last submesh, last nonNaN entry since it must monotonically increase
matrix_props.n_col = Mesh(end).indices.glob.unq_bounds.vert(2);

if strcmp(input.problemtype,'freeswim')
    switch input.bugtype
        case 'bacteria'
            switch input.tail.motorBC
                case 'freq'
                    matrix_props.n_rows = matrix_props.n_col * 3 + 6;  %usual unknowns plus 3 translation components and 3 rotation components of body
                case 'torque'
                    matrix_props.n_rows = matrix_props.n_col * 3 + 7;  %usual unknowns plus 3 translation components and 3 rotation components of body + rotation rate of tail
            end
        case {'dino' , 'sheet'}
            matrix_props.n_rows = matrix_props.n_col * 3 + 6;  %usual unknowns plus 3 translation components and 3 rotation components of body
    end
else
    matrix_props.n_rows = matrix_props.n_col * 3;
end
matrix_props.n_cols = matrix_props.n_rows; %number of columns in A matrix
matrix_props.Col_inds = save_Col_inds_mexed(Mesh,input.performance.nthreads);
