function [A_BIE_1,  A_BIE_2,  A_BIE_3,  A_BIE_4, RHS_BIE] = matrix_assembly_mex(Mesh,matrix_props,index_mapping,node_parameters,assembly_input)

%assembles all entries of A_BIE matrix block for both resistance problem and force-free cases

%if debug_mode is true, debug_info can be output

%produces up to 4 submatrices A_BIE_1 - 4 of full A_BIE to avoid mex error
%due to A_BIE being too big.  As a workaround, we insert A_BIE_1 - 4 into full A in
%matrix_assembly_mex_wrapper, in native Matlab, which has no dumb
%limitation on matrix size.

[refpoint, motor_orientation] = get_rotational_references(Mesh, assembly_input);
tail_ind = 0;
for i = 1:length(Mesh)
    if "Tail" == Mesh(i).name  % [Mesh.name] not allowed for code generation, hence this loop
        tail_ind = i; break;
    end
end

% temp = index_mapping;
% index_mapping = [];
% index_mapping.local_node2global_node = struct2cell(temp.local_node2global_node);
% index_mapping.global_node2local_node = reshape( struct2cell(temp.global_node2local_node) , [] ,1);
% index_mapping.global_node2global_u = temp.global_node2global_u;
% index_mapping.global_u2global_node = temp.global_u2global_node;

index_mapping = struct('local_node2global_node',{struct2cell(index_mapping.local_node2global_node)},...
    'global_node2local_node',{reshape( struct2cell(index_mapping.global_node2local_node) , [] ,1)},...
    'global_node2global_u',index_mapping.global_node2global_u,...
    'global_u2global_node',index_mapping.global_u2global_node);

% index_mapping.local_node2global_node = struct2cell(index_mapping.local_node2global_node);
% index_mapping.global_node2local_node = struct2cell(index_mapping.global_node2local_node)';


coder.extrinsic('tic');
coder.extrinsic('toc');
coder.extrinsic('num2str');

% debug_info.abs_error = [];
% debug_info.rel_error = [];
% debug_info.fun_evals = [];
% debug_info.flag = [];
% debug_info.min_value = [];
% debug_info.max_value = [];

% check to see if this is really needed for single submesh runs
% if length(Mesh) > 1 %multiple submeshes input, use global indices
%     tot_elems = Mesh(end).indices.glob.unq_bounds.elem(2);
% else %only one submesh input, use local indices
%     tot_elems = Mesh(1).n_elem;
% end


tic
%#codegen
%%
%preallocate A_BIE matrix and temporary submatrices needed for parallel assembly

%test anticipated size of A_BIE against dumb coder max array size
%limit and if problemo, split into sub-arrays
% splitting is based on columns to facilitate less parfor pain

%right now limited to a further division by 4 but this is easily expanded to
%whatever is needed

%numels_max = double(intmax('int32'));  %dumb Coder limit on numel of any array

numel_A_BIE = 3 * matrix_props.n_collocation * matrix_props.n_cols; % numel of matrix block corresponding to BIEs
% we get one vector equation (x y z components) by writing the BIE with x0 at a mesh node, and we do this for every node
% need to revise next comment, going back to x y z as 3rd dim due to parfor woes (loop indexing must be the same everywhere within loop)
% A_BIE will be organized with all x components, then all y, then all z (doing xyz xyz xyz ... is much harder to do with parfor)

% y_start = matrix_props.n_collocation; % row index of first y component of BIEs
% z_start = matrix_props.n_collocation*2;  %parfor doesn't allow putting these directly into indexing expressions...

assembly_input.performance.rigid_body_matrix_rotation = false;


persistent A_BIE0_1 A_BIE0_2 A_BIE0_3 A_BIE0_4
if assembly_input.performance.rigid_body_matrix_rotation
    if isempty(A_BIE0_1)
        firstrun = true;
        orientations0 = NaN(3,3,length(Mesh));
        for i = 1:length(Mesh)
            orientations0(:,:,i) = Mesh(i).orientation;
        end
    else
        firstrun = false;
    end
else
%     A_BIE0_1 = NaN; A_BIE0_2 = NaN; A_BIE0_3 = []; A_BIE0_4 = [];
    firstrun = true;
end


if numel_A_BIE < assembly_input.performance.numels_max  % no splits needed
    n_splits = 0;
    A_BIE_1 = zeros(matrix_props.n_collocation,3,matrix_props.n_cols); % x y z components of BIE as 2nd dim to faciliate parfor (which only allows the same
    % extremely restrictive (e.g. no : allowed) indexing expression for A_BIE_1 throughout the parfor loop, so no direct way to fill in x, y, z components as different rows)
    % weird to have just 3 columns and n_cols pages as 3rd dim but this avoids a presumably costly permute() to eventually make a 2D matrix (permute makes
    % a copy in memory, reshape does not)
    A_BIE_2 = [];    A_BIE_3 = [];     A_BIE_4 = [];
    %      A_BIE_2 = zeros(matrix_props.n_collocation,3,matrix_props.n_cols - size(A_BIE_1,2));
    if assembly_input.performance.rigid_body_matrix_rotation
        
        if firstrun
            switch assembly_input.performance.rigid_body_matrix_rotation
                case 'sparse'
                    A_BIE0_1_temp = sparse(matrix_props.n_collocation,3,matrix_props.n_cols); % BIE0 signifies the values at the initial orientation, to be
                    % rotated at later times/phases instead of recomputed from scratch
                    % a temp version is used for these because we need x, y, z as 3rd dim to initially fill in entries in the parfor but later need x, y, z
                    % as rows to do the tensor rotations and Coder doesn't allow reassignment of persistent variables (e.g. A_BIE0_1 = reshape(A_BIE0_1) )
                    % so we'll do initial creation with A_BIE0_temp and then assign A_BIE0 = reshape(A_BIE0_temp) and then use A_BIE0 in later function
                    % calls
                case 'full'
                    A_BIE0_1_temp = zeros(matrix_props.n_collocation,3,matrix_props.n_cols);
            end
            A_BIE0_2_temp = [];    A_BIE0_3_temp = [];     A_BIE0_4_temp = [];
        end
    end
    
elseif ceil(numel_A_BIE / 2)  < assembly_input.performance.numels_max  %split A_BIE into two parts
    n_splits = 1;
    A_BIE_1 = zeros(matrix_props.n_collocation,3, 3*round( ceil(matrix_props.n_cols / 2) / 3)); % split approx. in half but make sure # columns is divisible
    % by 3 so we don't split up an x y z triplet of components of traction or free-slip velocity
    A_BIE_2 = zeros(matrix_props.n_collocation,3,matrix_props.n_cols - size(A_BIE_1,2));
    A_BIE_3 = [];     A_BIE_4 = [];
    if assembly_input.performance.verbose
        disp('Dividing A_BIE into 2 parts as Coder limitation workaround.')
    end
    
    if assembly_input.performance.rigid_body_matrix_rotation
        
        if firstrun
            switch assembly_input.performance.rigid_body_matrix_rotation
                case 'sparse'
                    A_BIE0_1_temp = sparse(matrix_props.n_collocation,3,3*round( ceil(matrix_props.n_cols / 2) / 3)); % BIE0 signifies the values at the initial orientation, to be
                    % rotated at later times/phases instead of recomputed from scratch
                    A_BIE0_2_temp = sparse(matrix_props.n_collocation,3,matrix_props.n_cols - size(A_BIE0_1_temp,2));
                case 'full'
                    A_BIE0_1_temp = zeros(matrix_props.n_collocation,3,matrix_props.n_cols);
                    A_BIE0_2_temp = zeros(matrix_props.n_collocation,3,matrix_props.n_cols - size(A_BIE0_1_temp,2));
            end
            A_BIE0_3_temp = [];     A_BIE0_4_temp = [];
        end
    end
    
else %ceil(numel_A_BIE / 4)  < assembly_input.performance.numels_max  % split A_BIE into 4 parts
    n_splits = 2;
    A_BIE_1 = zeros(matrix_props.n_collocation,3,3*round(ceil(matrix_props.n_cols / 4)/3));
    A_BIE_2 = zeros(matrix_props.n_collocation,3,3*round(ceil(matrix_props.n_cols / 4)/3));
    A_BIE_3 = zeros(matrix_props.n_collocation,3,3*round(ceil(matrix_props.n_cols / 4)/3));
    A_BIE_4 = zeros(matrix_props.n_collocation ,3,matrix_props.n_cols - (size(A_BIE_1,2)+size(A_BIE_2,2)+size(A_BIE_3,2)));
    if assembly_input.performance.verbose
        disp('Dividing A_BIE into 4 parts as Coder limitation workaround.')
    end
    %code should fail if we end up needing further subdivisions -
    %can just add code for division into 8 parts etc...
    
    if assembly_input.performance.rigid_body_matrix_rotation
        
        if firstrun
            switch assembly_input.performance.rigid_body_matrix_rotation
                case 'sparse'
                    A_BIE0_1_temp = sparse(matrix_props.n_collocation,3,3*round(ceil(matrix_props.n_cols / 4)/3));
                    A_BIE0_2_temp = sparse(matrix_props.n_collocation,3,3*round(ceil(matrix_props.n_cols / 4)/3));
                    A_BIE0_3_temp = sparse(matrix_props.n_collocation,3,3*round(ceil(matrix_props.n_cols / 4)/3));
                    A_BIE0_4_temp = sparse( matrix_props.n_collocation ,3,matrix_props.n_cols - (size(A_BIE0_1_temp,2)+size(A_BIE0_2_temp,2)+size(A_BIE0_3_temp,2)));
                case 'full'
                    A_BIE0_1_temp = zeros(matrix_props.n_collocation,3,3*round(ceil(matrix_props.n_cols / 4)/3));
                    A_BIE0_2_temp = zeros(matrix_props.n_collocation,3,3*round(ceil(matrix_props.n_cols / 4)/3));
                    A_BIE0_3_temp = zeros(matrix_props.n_collocation,3,3*round(ceil(matrix_props.n_cols / 4)/3));
                    A_BIE0_4_temp = zeros( matrix_props.n_collocation ,3,matrix_props.n_cols - (size(A_BIE0_1_temp,2)+size(A_BIE0_2_temp,2)+size(A_BIE0_3_temp,2)));
            end
            
        end
    end
    
end

% if matrix_props.n_unknown_u > 0 % some verts have a free slip BC
% u_cols = zeros(matrix_props.n_col, matrix_props.n_unknown_u*3 ,3); % extra columns for coeffs corresponding to unknown u components in DL integrals
% u_rows = zeros(matrix_props.n_unknown_u*3 , matrix_props.n_cols);
% this is now part of BI_coll variable, lumped with traction

% end

RHS_BIE = zeros(matrix_props.n_collocation,3,matrix_props.n_RHS);
% need to revise for x y z organization% to match A_BIE, organized as all x, all y, all z, and 2nd dim is for possible multiple RHS

persistent RHS_BIE0
if assembly_input.performance.rigid_body_matrix_rotation
    if firstrun
        RHS_BIE0_temp = zeros(matrix_props.n_collocation,3,matrix_props.n_RHS);
    end
else
    RHS_BIE0_temp = [];
end

if assembly_input.performance.debug_mode
    %     abs_error_1 = NaN(matrix_props.n_col, tot_elems);
    %     rel_error_1 = NaN(matrix_props.n_col, tot_elems);
    %     fun_evals_1 = NaN(matrix_props.n_col, tot_elems);
    %     flag_1 =      NaN(matrix_props.n_col, tot_elems);
    %     min_value_1 = NaN(matrix_props.n_col, tot_elems);
    %     max_value_1 = NaN(matrix_props.n_col, tot_elems);
    %     time_1 =      NaN(matrix_props.n_col, tot_elems);
    
else
    %     abs_error_1 = [];
    %     rel_error_1 = [];
    %     fun_evals_1 = [];
    %     flag_1 =      [];
    %     min_value_1 = [];
    %     max_value_1 = [];
    %     time_1 = [];
end
if assembly_input.performance.verbose
    disp('initialized A_BIE');
end
% A_BIE = zeros(matrix_props.n_rows, matrix_props.n_cols); %includes extra variables and equations for mobility problem case
% disp('initialized A_BIE'); %assembling full A_BIE will have to be un-mexed for large dino meshes


coder.varsize('reference_nodes',[2 3 Inf],[true true true]); %for mex rules
reference_nodes = [0 1 0; 0 0 1];  %reference triangle normalized node coords

rule = 2;  %1st rule has ridiculously innaccurate (too conservative) error estimates, so never converges to the tolerances
[rule_constants.G, rule_constants.Weights, rule_constants.PTS] =  SMPRMS( 2, rule );  %always 2 dimensions and constant rule # for all integrals (2nd arg is rule #)

% [~, elem_range] = bounds_from_global_inds(Mesh);  %make it easier to find current submesh inside parfor
% [i_mesh_elems, local_elems] = global2local(1:tot_elems, Mesh, elem_range, 'elem');    %compute mesh_ind_elem and local_elem, indices to current element's submesh and local element index

% skipping traction integrals currently broken or more accurately, we could still recycle traction integrals easily but if we have free slip nodes, we
% can't easily deal with the additional U, Omega matrix entries going from resistance to mobility - would have to save int phiTn contributions from
% every node directly I guess?
% if ~assembly_input.skip_traction_integrals %if this is on, we skip the entire core part of A_BIE since we can copy it from the previous run (if we are doing mobility problem now and just did resistance problem)
%     Local_verts = cell(length(Mesh),1);  not used anywhere?


%     for mesh_ind_coll = 1:length(Mesh) %only parfor over one submesh at a time to avoid load balancing problem when skip_rigid_integrals = true

%         collocation_inds_submesh =  Mesh(mesh_ind_coll).indices.glob.unq_bounds.vert(1) : Mesh(mesh_ind_coll).indices.glob.unq_bounds.vert(2);  %global indices of collocation points for current submesh
%         [~, temp] = global2local(collocation_inds_submesh, Mesh, vert_range, 'vert');
%         local_verts = NaN(Mesh(mesh_ind_coll).indices.glob.unq_bounds.vert(2), 1);
%         local_verts(collocation_inds_submesh) = temp;  %slot in actual values into larger local_verts vector so we can directly index into it using col_i and avoid parfor errors
% %         Local_verts{mesh_ind_coll} = local_verts;  % not used anywhere?
%         mesh_loop_tic = tic;
%         eps2 = assembly_input.accuracy.eps2(mesh_ind_coll); % eps^2 going with reg. Stokeslets for all collocation points on current submesh

% was the below but write as col_inds unless parfor gets bitchy?
%         parfor (global_node_ind = Mesh(mesh_ind_coll).indices.glob.unq_bounds.vert(1) : Mesh(mesh_ind_coll).indices.glob.unq_bounds.vert(2)
%         parfor (global_node_ind = collocation_inds_submesh , assembly_input.performance.nthreads)  %rows of A_BIE, i.e. vertices

%     RHS_BIE0_temp = RHS_BIE0;
%     A_BIE0_1_temp = A_BIE0_1;  A_BIE0_2_temp = A_BIE0_2;  A_BIE0_3_temp = A_BIE0_3;  A_BIE0_4_temp = A_BIE0_4;


parfor (global_node_ind = 1:matrix_props.n_collocation , assembly_input.performance.nthreads)  %rows of A_BIE in sets of 3, i.e. x y z components of BIE for each node
    % this loops over collocation points
%         for global_node_ind = 1:matrix_props.n_collocation    %randperm(matrix_props.n_collocation)
    % instead of randomizing indices of actual meshes, why not introduce a vector that is a random permutation of 1:n_collocation and then here, use
    % randvec(global_node_ind) to get the actual collocation pt index that we get vert coords etc for.  This randomizes order of loop over collocation pts but
    % keeps meshes in original (presumably somewhat logical) node order.  After matrix solve, can un-randomize order of solution vector so that
    % traction, free slip velocity, etc match the original meshes.
    
    %                          for global_node_ind = Mesh(mesh_ind_coll).indices.glob.unq_bounds.vert(1) : Mesh(mesh_ind_coll).indices.glob.unq_bounds.vert(2)
    % col_i
    % for col_i = 1:matrix_props.n_col
    % col_i / matrix_props.n_col
    
    % global2local is a n_global_nodes cell array, each cell is a matrix with rows = submesh, local node ind, element ind, position in element
    %         index_mapping.global_node2local_node{global_ind}(end+1:end+size(element_inds,1),:) = [repmat([submesh_ind, local_node_ind],size(element_inds,1),1), element_inds, position_inds];
    
    temp = unique( index_mapping.global_node2local_node{global_node_ind}(:,[1 2]), 'rows' ); % submeshes and local node inds for this global node
    coll_local = struct('submeshes',temp(:,1),'nodes',temp(:,2));
    %     coll_local.submeshes = temp(:,1);
    %     coll_local.nodes = temp(:,2);
    
    
    %     for j = 1:length(Mesh)
    %         coll_i_local = find(global_node_ind == Mesh(j).indices.glob.vert,1);
    %         if ~isempty(coll_i_local)
    %             mesh_ind_coll = j; % index of submesh that current collocation pt belongs to (theoretically arbitrarily chosen between multiple submeshes
    %             % that share the coll pt)
    %             break
    %         end
    %         if isempty(coll_i_local) && j == length(Mesh)
    %             error('Can''t find global collocation index in any submesh');
    %         end
    %     end
    
    %     eps2 = assembly_input.accuracy.eps2(mesh_ind_coll); % eps^2 going with reg. Stokeslets for all collocation points on current submesh
    % innertic = tic;
    
    AE = [];   NV = [];     FL = logical([]); %stop parfor complaints about temp variables
    r_nodes = []; u_element_nodes = []; RHS_coll0 = []; BI_coll0 = [];
    %             local_vert = local_verts(global_node_ind);  replaced by below, if it doesn't work may have to go back to original
    %             [~, coll_i_local] = global2local(global_node_ind, Mesh, vert_range, 'vert');
    
    
    u_col = reshape(  squeeze(Mesh(coll_local.submeshes(1)).u(coll_local.nodes(1),:,:)), 3 , []  ); % needed if including DL  3 x flowcases
    % reshape needed in case there's only one flowcase, to make sure 
    % arbitrarily look up u_col for first submesh/local node this global node appears in - u_col should be the same for all appearances
    
    [full_DL_calcs, inds] = setdiff( coll_local.submeshes , find(assembly_input.performance.eliminate_DL | assembly_input.performance.DL_singularity_removal == 1) ); % submeshes that this collocation pt is on for which we are not eliminating the DL or using trick
    % note, eliminate_DL should always be either true or false while DL_singularity_removal can be NaN for bodies that should always have DL eliminated
    % e.g. sheets
    alphas = NaN(length(full_DL_calcs),1);
    for n = 1:length(full_DL_calcs)
        alphas(n) = Mesh( full_DL_calcs(n) ).solid_angle(coll_local.nodes(inds(n)));
    end
    
    alpha_factor = sum(alphas) + (1 - length(full_DL_calcs))*4*pi; % this is what multiplies the collocation velocity u(x_0) in the BIE; in the typical case
    % of u(x_0) being a member of just one body, alpha_factor = alpha.  For e.g. u(x_0) on dino body and transverse, since transverse alpha = 4*pi, we
    % still get alpha of just the body here.
    
    %     if assembly_input.performance.eliminate_DL(coll_local.submesh) || assembly_input.performance.DL_singularity_removal
    %         alpha = 4*pi; % whether this mesh is a sheet, a rigid body, or a deforming body (that we don't need correct traction for), eliminating the
    %         %DL results in solid angle effectively becoming 4*pi regardless of the actual geometric value
    %         %keeping the DL integral but using the singularity removal identity for the mesh containing the col pt also results in alpha effectively becoming 4*pi
    %     else
    %         alpha = Mesh(coll_local.submesh).solid_angle(coll_local.node);  % solid angle of surface at current collocation pt:  2*pi for smooth, closer to zero for a valley, closer to 4*pi for a ridge or peak
    %     end
    
    %  [DL_elem, ~, ~, ~] = adsimp( 2, reference_nodes, 54,  assembly_input.accuracy.integration_tol.stresslet.phiTn.maxevals, abstol.stresslet.phiTn,...
    %      assembly_input.accuracy.integration_tol.stresslet.phiTn.reltol ,rule  ,...
    %      rule_constants ,element_nodes,shape_parameters,integrand_constants,'reg_stresslet',{'phiTn'},assembly_input.performance.DL_singularity_removal);
    %
    
          x_col = Mesh(coll_local.submeshes(1)).nodes(coll_local.nodes(1),:)';      
                
    integrand_constants = struct('eps2',NaN,...
        'x_col',x_col,...
        'element_nodes',NaN(6,3),...
        'shape_parameters',NaN(3,1),...
        'integral_type',"",...
        'DL_singularity_removal',false,... % doesn't matter for single layer integral, will overwrite later for double layer integral
        'refpoint',NaN(3,1),...
        'field_vel',NaN(6,3)); % either 1 - 6 or 0 for if/where collocation pt is on current element, only matters for DL singularity removal trick with DL integrals
    % as with u_col, shouldn't matter which instance of this node we look up the coordinates for, should all be the same
    % eps is allowed to vary between submeshes so it is defined later in the loop over elements
    % 'coll_node_ind', [],...
    
%     integrand_constants = struct;
%     integrand_constants.x_col = x_col;
%     integrand_constants.DL_singularity_removal = false;  % doesn't matter for single layer integral, will overwrite later for double layer integral
    
    r_col = x_col - refpoint; % vector from refpoint to collocation pt
    %     BC_type_col = Mesh(coll_local.submesh).vert_BC_type(coll_local.node); % may be needed if including DL
    BC_type_col = node_parameters.BC_type(global_node_ind);
    %     global_u_ind = Mesh(coll_local.submesh).indices.glob.unknown_u(coll_local.node); % index for the unknown velocity at the current collocation point,
    % between 1 - n where n is number of free slip nodes
    global_u_ind = index_mapping.global_node2global_u( global_node_ind );
    
    
    
    %%temporary matrix for boundry integrals (BI_coll), row is x, y, z components for current collocation point BIE and columns are coeffs for traction components at each vert, to later insert into submatrices
    BI_coll = zeros(3,matrix_props.n_cols); % 3 rows for x, y, z components of BIE, columns for x y z components of traction followed by x y z components of unknown u (for free slip coll pts) followed by any unknown kinematic variables (U, Omega, omega)
    RHS_coll = zeros(3,matrix_props.n_RHS); % just for current collocation pt, for all flowcases
    
    if assembly_input.performance.rigid_body_matrix_rotation && firstrun
        BI_coll0 = zeros(3,matrix_props.n_cols);
        RHS_coll0 = zeros(3,matrix_props.n_RHS);
    end
    
    % there may be a small penalty to having only 3 rows and many columns (since arrays are stored column-major) but seems unlikely this is significant
    % since these arrays aren't *that* big?  Can check where time is spent with Profiler before mexing?
    % also, if we flipped the dimensions for these, we'd have to do a transpose later to insert them into the larger arrays, which would be
    % expensive, so it's a matter of which penalty is larger...
    
    %A_y_temp = zeros(matrix_props.n_col*3,1); A_z_temp = zeros(matrix_props.n_col*3,1);
    %             A_temp = BI_coll; %what actually gets inserted into A - is the same as BI_coll for no-slip BCs but will be different for free-slip
    
    %for debugging + mex rules
    %     abs_error_temp = [];
    %     rel_error_temp = [];
    %     fun_evals_temp = [];
    %     flag_temp = [];
    %     min_value_temp = [];
    %     max_value_temp = [];
    %             time_temp = [];
    
    if assembly_input.performance.debug_mode
        %                 integral = NaN(3,18,tot_elems);  abserror = integral;
        %                 numevals = NaN(1,tot_elems);
        %                 exitflag = numevals;
        
        %         abs_error_temp = NaN(1,tot_elems);
        %         rel_error_temp = NaN(1,tot_elems);
        %         fun_evals_temp = NaN(1,tot_elems);
        %         flag_temp =      NaN(1,tot_elems);
        %         min_value_temp = NaN(1,tot_elems);
        %         max_value_temp = NaN(1,tot_elems);
        %                 time_temp = NaN(1,tot_elems);
        
    end
    
    %             elem_t1 = [];  %appease parfor rules
    
    
    
    %         RHS_temp0 = zeros(3,matrix_props.n_RHS);
    for mesh_ind_elem = 1:length(Mesh)
        % ignore interaction if the element submesh is not a member of any submeshes the coll pt is on
        if assembly_input.accuracy.ignore_interaction && ~ismember( mesh_ind_elem, coll_local.submeshes )
            
            %collocation point and element are in different submeshes - skip integrals
            continue
        end
        
        % if coll pt is a member of the element-containing submesh, then can do tensor rotations since the coll pt and element must be rigidly connected, even if
        % another submesh the coll pt belongs to is deforming
        if assembly_input.performance.rigid_body_matrix_rotation && ~firstrun && Mesh(mesh_ind_elem).is_mesh_rigid && ismember( mesh_ind_elem , coll_local.submeshes )
            % we can simply rotate the tensor contributions coming from the nodes on this element instead of integrating again
            continue
        end
        
        integrand_constants.eps2 = assembly_input.accuracy.eps2(mesh_ind_elem); % each submesh has an eps corresponding to all its elements
        % since we're assuming eps is piecewise constant across elements, there can be discontinuities between submeshes, but this should be OK?
        
        for local_elem = 1:Mesh(mesh_ind_elem).n_elements
            %     for elem_i = 1:tot_elems  %elem_i is a global index
            %         elem_i / Mesh.n_elem
            %                 if assembly_input.performance.debug_mode
            %                     elem_t1 = clock;
            %                 end
            
            %                 [mesh_ind_elem, local_elem] = global2local(elem_i, Mesh, elem_range, 'elem');    %compute mesh_ind_elem and local_elem, indices to current element's submesh and local element index
            %         mesh_ind_elem = i_mesh_elems(elem_i);   local_elem = local_elems(elem_i);
            
            %          ind2 = matrix_props.global2local.elem.permuted_inds(elem_i);  % probably a dummy operation since elem_i = ind2 unless you decide to randomly permute elements as well as verts
            %     mesh_ind_elem = matrix_props.global2local.elem.submesh(ind2); %submesh index
            %     local_elem = matrix_props.global2local.elem.local(ind2);  %element index
            
            
            BI_elem = zeros(3,matrix_props.n_cols); % stores cumulative matrix contributions for integrals over a single element
            RHS_elem = zeros(3,matrix_props.n_RHS); % stores cumulative RHS contributions for integrals over a single element
            
            
            %want mesh_ind_elem below since we are generally interested in current
            %element; node only matters as far as it's distance away when computing reg stokeslet
            %function
            element_node_inds = Mesh(mesh_ind_elem).elements(local_elem,:); % local inds for nodes of current element  1 x 6
            element_nodes = Mesh(mesh_ind_elem).nodes(element_node_inds,:); % coords of nodes of current element   6 x 3
            element_global_node_inds = index_mapping.local_node2global_node{mesh_ind_elem}(element_node_inds); % global inds for nodes of current element  6 x 1
            
            %         element_global_node_inds = Mesh(mesh_ind_elem).indices.glob.vert(element_node_inds); % global inds for verts of current element
            %         submesh_members = {Mesh(mesh_ind_elem).indices.submesh_members{element_node_inds}}; % submeshes each vert of this element belong to
            
            
            BC_types = node_parameters.BC_type( element_global_node_inds ); % BC for each element node:  1 for no slip, 2 for free slip
            if ~assembly_input.performance.eliminate_DL(mesh_ind_elem) % include DL integrals over current submesh elements?
                u_element_nodes = ( Mesh(mesh_ind_elem).u(element_node_inds,:,:) ); % known u contributions (body frame velocities for mobility problem, fixed frame velocities
                % resistance problem) at each element node.   6 x 3 x F flowcases
            end
            
            if strcmp( assembly_input.problemtype , 'mobility' )
                r_nodes = element_nodes' - repmat(refpoint,1,6); % 3 x 6 vectors from refpoint to each vert of current element
            end
            
            unknown_u_global_inds = index_mapping.global_node2global_u(element_global_node_inds); % indices into matrix in case there are unknown fluid u due
            % to free slip BC at some nodes of current element
            
            integrand_constants.element_nodes = element_nodes;
            integrand_constants.shape_parameters = Mesh(mesh_ind_elem).shape_parameters(local_elem,:)';  %[alpha beta gamma]'  3 x 1
             
            integrand_constants.coll_node_ind = find( global_node_ind == element_global_node_inds);  % needed to appease Coder, which requires this be scalar
%                 if ~isempty(temp)
%                     temp = temp(1);
%                 end
%                 integrand_constants.coll_node_ind = temp;  % *really* local index (1-6) of node of current element matching current collocation pt, otherwise empty
               
            %         col_inds = matrix_props.Col_inds(elem_i,:); %Col_inds was defined based on global indices so use elem_i not local_elem
            % these are the column inds of A_BIE that go with the traction (x y z components) at each of 6 nodes of the current element
            
            
            %scale abs integration tolerance by element area, since big elements
            %should be allowed more total error than small elements
            abstol_stokeslet =  assembly_input.accuracy.integration_tol.stokeslet.abstol * Mesh(mesh_ind_elem).area(local_elem);
            %don't scale reltol, since reltol should already account for area
            %since the integral itself should scale with area
            
            %                 abstol.stresslet.Tn =  assembly_input.accuracy.integration_tol.stresslet.Tn.abstol * Mesh(mesh_ind_elem).area(local_elem);
            %                 abstol.stresslet.rTn =  assembly_input.accuracy.integration_tol.stresslet.rTn.abstol * Mesh(mesh_ind_elem).area(local_elem);
            
            %the important part:  adaptively compute boundary integrals of reg stokeslet function over
            %current element, with respect to current collocation point
            integrand_constants.integral_type = "reg_stokeslet";
            
            if assembly_input.performance.debug_mode
                % shatlab hangs at almost no CPU usage with the tic/toc
                % inside adsimp so timing inside the parfor doesn't work even though it
                % compiles without errors
                %[SL_elem, AE, NV, FL, Time] = adsimp( 2, reference_nodes, 54,  assembly_input.accuracy.integration_tol.traction.maxevals,abstol, assembly_input.accuracy.integration_tol.traction.reltol ,rule  , rule_constants ,element_nodes,shape_parameters,integrand_constants,'reg_stokeslet', true);
                [SL_elem, AE, NV, FL] = adsimp( 2, reference_nodes, 54,  assembly_input.accuracy.integration_tol.stokeslet.maxevals,abstol_stokeslet, assembly_input.accuracy.integration_tol.stokeslet.reltol ,rule  , rule_constants ,integrand_constants);
                %                     SL_elem = ones(54,1);  AE = ones(54,1);  NV = 1;  FL = 1;
            else %don't bother storing extra outputs
                %                     [SL_elem, ~, ~, ~, ~] = adsimp( 2, reference_nodes, 54,  assembly_input.accuracy.integration_tol.traction.maxevals,abstol, assembly_input.accuracy.integration_tol.traction.reltol ,rule  , rule_constants ,element_nodes,shape_parameters,integrand_constants,'reg_stokeslet', false);
                
                [SL_elem, ~, ~, ~] = adsimp( 2, reference_nodes, 54,  assembly_input.accuracy.integration_tol.stokeslet.maxevals,abstol_stokeslet, assembly_input.accuracy.integration_tol.stokeslet.reltol ,rule  , rule_constants ,integrand_constants);
                %                 SL_elem = ones(54,1);  AE = ones(54,1);  NV = 1;  FL = 1;
                SL_elem = SL_elem * 1/2 / assembly_input.constants.mu;  % assumes that solid angle (2*pi for a smooth region) multiplies u_c on left side of BIE
                
            end
            %               out =  hS*[S(1,1)*phi; S(1,2)*phi; S(1,3)*phi; ...
            %             S(2,1)*phi; S(2,2)*phi; S(2,3)*phi; ...
            %             S(3,1)*phi; S(3,2)*phi; S(3,3)*phi];
            % integrals of phiS
            submat = reshape(SL_elem,18,3)'; % 6 nodes, 3 x y z components of traction per node = 3 rows for x y z components of BIE,
            % 18 columns for x, y, z traction components and 6 phi
            
            %entries in matrix are total boundary integrals (BI) so can do a cumulative sum
            inds = ( repmat( 3 * (element_global_node_inds - 1),1,3) + repmat([1 2 3],6,1) ); % indices of matrix columns corresponding to x,y,z traction components at global nodes of
            % current element   6 x 3 each column are the matrix column inds of the 6 nodes, and we have 3 columns for x, y, z components of traction
            
            
            BI_elem(:,inds(:)) = BI_elem(:,inds(:)) + submat;
            %                 BI_elem(2,col_inds) = BI_elem(2,col_inds) + submat(:,2)';
            %                 BI_elem(3,col_inds) = BI_elem(3,col_inds) + submat(:,3)';
            
            %             if assembly_input.performance.rigid_body_matrix_rotation && firstrun && Mesh(mesh_ind_elem).rigid && mesh_ind_elem == mesh_ind_coll
            %                 BI0(:,col_inds) = BI0(:,col_inds) + submat';  % this version of BI_coll will only have constributions from integrals for elements and coll
            %                 % pts on the same rigid body
            %             end
            
            
            if ~assembly_input.performance.eliminate_DL(mesh_ind_elem)  % compute DL integrals over this element
                integrand_constants.integral_type = "reg_stresslet";
                integrand_constants.DL_singularity_removal = assembly_input.performance.DL_singularity_removal(mesh_ind_elem);
                abstol_stresslet =  assembly_input.accuracy.integration_tol.stresslet.phiTn.abstol * Mesh(mesh_ind_elem).area(local_elem);
                % integrating phiTn regardless of no slip or free slip, but need to handle all flowcases for no slip
                [DL_elem, ~, ~, ~] = adsimp( 2, reference_nodes, 54,  assembly_input.accuracy.integration_tol.stresslet.phiTn.maxevals, abstol_stresslet, assembly_input.accuracy.integration_tol.stresslet.phiTn.reltol ,rule  , rule_constants , integrand_constants);
                DL_elem = DL_elem * 1/2;  %  assumes that solid angle (2*pi for a smooth region) multiplies u_c on left side of BIE
                
                switch assembly_input.problemtype
                    case "resistance"
                        
                        %                                 assembly_input.BC_type(mesh_ind_elem) % only matters for inserting values into A_BIE or RHS, either way we're integrating phiTn
                        % mobility:  for elements that are shared between no slip, free slip vertices, can simply compute all variations of the integrand (e.g.
                        % Tn, rTn, phiTn) and just only use the ones needed for each vertex of the element
                        
                        
                        
                        % each Mesh(i).u is Nx3xM where M is length(flowcases), N is #verts
                        %                             u_cols = zeros(matrix_props.n_col, matrix_props.n_unknown_u*3 ,3); % extra columns for coeffs corresponding to unknown u components in DL integrals
                        %                             u_rows = zeros(matrix_props.n_unknown_u*3 , matrix_props.n_cols);
                        
                        
                        
                        for vi = 1:6
                            % phiTn for phi of current vert
                            phiTn = [ DL_elem(vi)  DL_elem(6+vi) DL_elem(12+vi) ;...
                                DL_elem(18+vi)  DL_elem(24+vi)  DL_elem(30+vi) ;...
                                DL_elem(36+vi) DL_elem(42+vi)  DL_elem(48+vi) ;  ];
                            
                            switch BC_types(vi)
                                
                                case 1 % no slip
                                    
                                    % flowcases, e.g. forced x y z translation and x y z rotation to calculate diffusivities
                                    if assembly_input.performance.DL_singularity_removal(mesh_ind_elem) && BC_type_col == 1
                                        % no slip, u_col is prescribed and known
                                        vel_mod = reshape( squeeze(u_element_nodes(vi,:,:)),3,[]) - u_col;  %u_element_nodes is 6 x 3 x flowcases, u_col is 3 x flowcases
                                        % free slip, u_col is unknown
                                    else % either not doing trick, or doing trick but u_col is unknown
                                        vel_mod = reshape( squeeze(u_element_nodes(vi,:,:)),3,[]);
                                    end
                                    % velmod is 3 x flowcases
                                    % in any case, add integral corresponding to known u_v to RHS
                                    RHS_elem = RHS_elem - phiTn * vel_mod; % RHS_elem is 3 x flowcases
                                    % minus since we're moving this term from LHS to RHS
                                    %                                     if assembly_input.performance.rigid_body_matrix_rotation && firstrun && Mesh(mesh_ind_elem).rigid && mesh_ind_elem == mesh_ind_coll
                                    %                                         RHS_temp0 = RHS_temp0 - phiTn * vel_mod;   % this version of BI_coll will only have constributions from integrals for elements and coll
                                    %                                         % pts on the same rigid body
                                    %                                     end
                                    
                                case 2 % free slip
                                    % in case of no trick, these coeffs go with u for this vert (across columns for each u component, down rows for each BIE
                                    % component)
                                    inds = matrix_props.n_collocation*3 + 3 * (unknown_u_global_inds(vi) - 1) + [1 2 3]; % corresponding to inds for 3 components of unknown u
                                    % matrix columns are organized with unknown u following unknown traction
                                    BI_elem(:,inds ) = BI_elem(:,inds ) + phiTn;
                                    
                                    % in case of trick, these coeffs still go with u for this vert, AND -coeffs go with u_col:
                                    if assembly_input.performance.DL_singularity_removal(mesh_ind_elem) && BC_type_col == 1
                                        % if u_col is no-slip, then add sum(coeffs .* u_col , 2) to RHS.
                                        % flowcases, e.g. forced x y z translation and x y z rotation to calculate diffusivities
                                        RHS_elem = RHS_elem + phiTn * u_col; % minus sign from moving from LHS to RHS negated since we want to add integral of -u_col
                                        
                                    end
                                    
                            end % no slip vs free slip
                            
                            if assembly_input.performance.DL_singularity_removal(mesh_ind_elem) && BC_type_col == 2
                                % if u_col is free-slip, also add integral corresponding to unknown u_c to coeffs in matrix
                                inds = matrix_props.n_collocation*3 + 3 * (global_u_ind - 1) + [1 2 3]; % corresponding to inds for 3 components of unknown u
                                % matrix columns are organized with unknown u following unknown traction
                                % minus sign since we were integrating (u_v - u_c)*Tn
                                BI_elem(:,inds ) = BI_elem(:,inds) - phiTn;
                            end
                            
                        end % element_nodes of current element
                        
                        
                        
                    case "mobility"
                        
                        
                        %mobility, noslip   same matrix size    add DL terms into coeffs for U, Omega, and export DL terms for u_r to RHS
                        % integrate Tn for U coeffs (9 values for 3 U components, 3 BIE components), rTn for Omega coeffs (9 values), phiTn to overall integrate uTn
                        % where u = known, interped u_r on flagella (6 x 9 = 54 values, combined with known u_r on element verts to yield 3 values for 3 components of BIE, exported to RHS)
                        % in trick, Tn is set to zeros since int(U - U_c)Tn will be zero everywhere
                        % in trick, we overall int( Omega x (r - r_c) Tn ) so that integrand goes to zero at collocation vertex
                        % in trick, we manually zero out phiTn when an element vertex is the same as the current collocation pt
                        %                                 maxevals_temp = [ repmat(assembly_input.accuracy.integration_tol.stresslet.Tn.maxevals,9,1); repmat(assembly_input.accuracy.integration_tol.stresslet.rTn.maxevals,9,1); repmat(assembly_input.accuracy.integration_tol.stresslet.phiTn.maxevals,54,1); ];
                        %                                 abstol_temp = [ repmat(abstol.stresslet.Tn,9,1); repmat(abstol.stresslet.rTn,9,1);  repmat(abstol.stresslet.phiTn,54,1); ];
                        %                                 reltol_temp = [ repmat(assembly_input.accuracy.integration_tol.stresslet.Tn.reltol,9,1); repmat(assembly_input.accuracy.integration_tol.stresslet.rTn.reltol,9,1); repmat(assembly_input.accuracy.integration_tol.stresslet.phiTn.reltol,54,1); ];
                        
                        %                                 [DL_elem, ~, ~, ~] = adsimp( 2, reference_nodes, 54,  maxevals_temp, abstol_temp,reltol_temp ,rule  , rule_constants ,element_nodes,shape_parameters,integrand_constants,'reg_stresslet',{'Tn','rTn','phiTn'},assembly_input.performance.DL_singularity_removal);
                        %                                 Tn = reshape(DL_elem(1:9),3,3);  rTn = reshape(DL_elem(10:18),3,3);  phiTn = DL_elem(19:72);
                        
                        for vi = 1:6
                            % phiTn for phi of current vert
                            phiTn = [ DL_elem(vi)  DL_elem(6+vi) DL_elem(12+vi) ;...
                                DL_elem(18+vi)  DL_elem(24+vi)  DL_elem(30+vi) ;...
                                DL_elem(36+vi) DL_elem(42+vi)  DL_elem(48+vi) ;  ];
                            
                            switch BC_types(vi)
                                
                                case 1 % no slip
                                    
                                    % keeping the ability to handle multiple "flowcases" here for giggles - i.e. if one wanted to try
                                    % relative surface velocities, perhaps for a squirmer? (geometry would need to be the same over time
                                    % for each case, which is kind of limiting....)
                                    if assembly_input.performance.DL_singularity_removal(mesh_ind_elem) && BC_type_col == 1
                                        % no slip, u_col is prescribed and known
                                        vel_mod = reshape( squeeze(u_element_nodes(vi,:,:)),3,[]) - u_col;  %u_element_nodes is 6 x 3 x flowcases, u_col is 3 x flowcases
                                        r_mod = r_nodes(:,vi) - r_col;
                                        U_factor = 0; % int(U_v - U_c)Tn is always zero, even when current node is not same as current coll pt, since U is
                                        % a constant everywhere
                                        % free slip, u_col is unknown
                                    else % either not doing trick, or doing trick but u_col is unknown
                                        vel_mod = reshape( squeeze(u_element_nodes(vi,:,:)),3,[]);
                                        r_mod = r_nodes(:,vi);
                                        U_factor = 1;
                                    end
                                    % vel_mod is 3 x flowcases
                                    
                                    % rphiTn are the coeffs going with (Omega x r)phiTn   3 x 3
                                    rphiTn = [r_mod(2)*phiTn(:,3) - r_mod(3)*phiTn(:,2) , r_mod(3)*phiTn(:,1) - r_mod(1)*phiTn(:,3) , r_mod(1)*phiTn(:,2) - r_mod(2)*phiTn(:,1) ];
                                    % may have to repmat each r into 3 rows for Coder
                                    
                                    % these entries hold whether or not we do the trick since the trick is taken care of when we compute Tn, rTn and with
                                    % U_factor and r_mod
                                    inds = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + (1:6);
                                    BI_elem(:,inds ) = BI_elem(:,inds ) + [repmat(U_factor,3,3).*phiTn  rphiTn]; % coeffs for U and Omega
                                    % need U_factor here to zero out the U phiTn term with trick even when current node is not same as current coll pt, since U phiTn
                                    % term always cancels out when doing trick
                                    
                                    if assembly_input.rotating_flagellum && strcmp(assembly_input.Tail.motorBC,'torque') && mesh_ind_elem == tail_ind
                                        % (only accounting for one tail for now, need to generalize for multiflagellated)
                                        % if element node is on bacterial tail and no slip, with unknown omega
                                        ind = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 7; % coeff for unknown omega for tail rotation relative to body
                                        BI_elem(:,ind ) = BI_elem(:,ind ) - rphiTn * motor_orientation; % see "double layer" work
                                        % negative sign due to how we define the contribution of relative tail rotation (see Schuech et al PNAS SI)
                                    end
                                    
                                    
                                    % in any case, add integral corresponding to known u_v to RHS
                                    RHS_elem = RHS_elem - phiTn * vel_mod; % RHS_elem is 3 x flowcases
                                    % minus since we're moving this term from LHS to RHS
                                    
                                case 2 % free slip
                                    % in case of no trick, these coeffs go with u for this vert (across columns for each u component, down rows for each BIE
                                    % component)
                                    inds = matrix_props.n_collocation*3 + 3 * (unknown_u_global_inds(vi) - 1) + [1 2 3]; % corresponding to inds for 3 components of unknown u
                                    % matrix columns are organized with unknown u following unknown traction
                                    BI_elem(:,inds ) = BI_elem(:,inds ) + phiTn;
                                    
                                    % in case of trick, these coeffs still go with u for this vert as included directly above, AND -coeffs also go with u_col:
                                    if assembly_input.performance.DL_singularity_removal(mesh_ind_elem) && BC_type_col == 1
                                        % need to include -int(U phiTn) - int(Omega x r phiTn) for collocation pt
                                        % rphiTn are the coeffs going with (Omega x r)phiTn   3 x 3
                                        rphiTn = [r_col(2)*phiTn(:,3) - r_col(3)*phiTn(:,2) , r_col(3)*phiTn(:,1) - r_col(1)*phiTn(:,3) , r_col(1)*phiTn(:,2) - r_col(2)*phiTn(:,1) ];
                                        % may have to repmat each r into 3 rows for Coder
                                        
                                        inds = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + (1:6);
                                        BI_elem(:,inds ) = BI_elem(:,inds ) - [phiTn  rphiTn]; % coeffs for U and Omega
                                        
                                        if assembly_input.rotating_flagellum && assembly_input.Tail.motorBC == "torque" && ismember( tail_ind , coll_local.submeshes )
                                            ind = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 7; % coeff for unknown omega for tail rotation relative to body
                                            BI_elem(:,ind ) = BI_elem(:,ind ) + rphiTn * motor_orientation; % see "double layer" work
                                            % plus:  normally would be minus due to how tail rotation is accounted for, see Schuech PNAS 2019, but we are
                                            % subtracting u_col via the trick so it becomes plus
                                        end
                                        
                                        
                                        % also add int(u_r phiTn) at collocation pt to RHS.
                                        % flowcases, e.g. forced x y z translation and x y z rotation to calculate diffusivities
                                        RHS_elem = RHS_elem + phiTn * u_col; % minus sign from moving from LHS to RHS negated since we want to add integral of -u_col
                                        
                                    end
                                    
                            end % no slip vs free slip
                            
                            if assembly_input.performance.DL_singularity_removal(mesh_ind_elem) && BC_type_col == 2
                                % if u_col is free-slip, also add integral corresponding to unknown u_c to coeffs in matrix
                                inds = matrix_props.n_collocation*3 + 3 * (global_u_ind - 1) + [1 2 3]; % corresponding to inds for 3 components of unknown u
                                % matrix columns are organized with unknown u following unknown traction
                                % minus sign since we were integrating (u_v - u_c)*Tn
                                BI_elem(:,inds ) = BI_elem(:,inds) - phiTn;
                            end
                            
                            
                        end % element_nodes of current element
                        
                        
                        
                end % resistance problem vs mobility problem
                
                
                % resistance problem, noslip     same matrix size   DL should be zero but can incl as RHS only, need to calc here and export out for every flowcase
                % integrate phiTn where overall aim is to integate uTn (no trick) or (u - u_c)Tn (trick) where u is known, interped u at any point, for all flowcases
                % in trick, we manually zero out phiTn when an element vertex is the same as the current collocation pt
                % 3 values for 3 components of BIE (for each flowcase)
                
                
                % resistance problem, freeslip  bigger matrix    DL nonzero, need to add all DL terms involving u as unknowns.  RHS = zero for BIEs and RHS = imposed (U+Omegaxr).n for additional equations
                % integrate phiTn for coeffs that multiply 3 components of 6 u_v velocities at verts, or in case of trick, (u_v - u_c) at each vert
                % 54 values from 6 phi * 9 Tn entries
                % leave u_c as an unknown in BIE, then add additional equations for u_c.n = known values
                
                
                %mobility problem, noslip   same matrix size    add DL terms into coeffs for U, Omega, and export DL terms for u_r to RHS
                % integrate Tn for U coeffs (9 values for 3 U components, 3 BIE components), rTn for Omega coeffs (9 values), phiTn to overall integrate uTn
                % where u = known, interped u_r on flagella (6 x 9 = 54 values, combined with known u_r on element verts to yield 3 values for 3 components of BIE, exported to RHS)
                % in trick, Tn is set to zeros since int(U - U_c)Tn will be zero everywhere
                % in trick, we overall int( Omega x (r - r_c) Tn ) so that integrand goes to zero at collocation vertex
                % in trick, we manually zero out phiTn when an element vertex is the same as the current collocation pt
                
                
                %mobility problem, freeslip  bigger matrix
                % integrate phiTn for coeffs that multiply 3 components of 6 u_v velocities at verts, or in case of trick, (u_v - u_c) at each vert
                % 54 values from 6 phi * 9 Tn entries
                % same integration as resistance problem, freeslip since we don't know any u_v velocities here either
                % difference from resistance problem comes in during matrix assembly, where additional constraint eqs on u_c.n must be expanded to
                % U.n + (Omega x r).n + u_r.n and we have the additional eqs for sum(forces), sum(torques) = 0 to solve for U, Omega
                
            end % include DL
            
            
            
            
            
            if assembly_input.performance.debug_mode
                
                %             abs_error_temp(elem_i) = max(abs(AE));  %the largest abs error over all 54 component integrals
                %             rel_error_temp(elem_i) = max(    abs(AE ./ SL_elem)  );  %the largest rel error over all 54 component integrals
                %             fun_evals_temp(elem_i) = NV;
                %             flag_temp(elem_i) = FL;
                %
                %             min_value_temp(elem_i) = min(abs(SL_elem) );  %smallest (absolute) component integral
                %             max_value_temp(elem_i) = max(abs(SL_elem) );  %largest (absolute) component integral
                
                %                     time_temp(elem_i) = Time;
                
            end
            
            %           BI_elems = BI_elems + BI_elem;
            %              RHS_elems = RHS_elems + RHS_elem;
            
            BI_coll = BI_coll + BI_elem;
            RHS_coll = RHS_coll + RHS_elem;
            
            if assembly_input.performance.rigid_body_matrix_rotation && firstrun && Mesh(mesh_ind_elem).is_mesh_rigid && ismember(mesh_ind_elem , coll_local.submeshes)
                
                BI_coll0 = BI_coll0 + BI_elem;
                RHS_coll0 = RHS_coll0 + RHS_elem;   % this version of BI_coll will only have constributions from integrals for elements and coll
                % pts on the same rigid body
            end
            
            
            
        end  % elements
        
    end % submeshes
    
    
    % currently redoing below contributions even for rigid body coll pts + elements since I'm not sure if/how the rotations would work for these.
    % the cost of redoing them is probably very small though
    
    % add matrix and/or RHS contributions from the alpha_factor*u_c = side of the BIEs, whether alpha_factor is typically the actual solid angle at the coll pt or set to
    % 4 pi due to eliminating the DL or using the singularity removal identity on the mesh containing the coll pt
    % note that if using singularity removal identity on any meshes NOT containing the coll pt, then the last term in the identity eq is zero for each of these DL integrals
    % (it is (alpha - 4*pi)u_c and that alpha is 4 pi since the coll point is not in/on that mesh - note, not same alpha_factor as the LHS of the BIE) so we
    % don't need to do anything about it here
    
    switch BC_type_col
        case 1 % no slip
            % RHS_elem is 3 x flowcases
            RHS_coll = RHS_coll + alpha_factor*u_col; % contribution of known u_r    u_col is 3 x flowcases
            
            switch assembly_input.problemtype
                
                %                         case 'resistance' % known u_c goes on RHS
                
                case "mobility" % u_c = U + (Omega x r) + u_r
                    inds = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + (1:6);
                    % subtract to move unknowns into matrix side of BIE
                    BI_coll(1,inds ) = BI_coll(1,inds ) - alpha_factor*[1 0 0   0          r_col(3)   -r_col(2)]; % coeffs from U and (Omega x r)
                    BI_coll(2,inds ) = BI_coll(2,inds ) - alpha_factor*[0 1 0  -r_col(3)   0           r_col(1)]; % coeffs from U and (Omega x r)
                    BI_coll(3,inds ) = BI_coll(3,inds ) - alpha_factor*[0 0 1   r_col(2)  -r_col(1)    0]; % coeffs from U and (Omega x r)
                    
                    if assembly_input.rotating_flagellum && assembly_input.Tail.motorBC == "torque" && ismember( tail_ind , coll_local.submeshes)
                        % (tail hardcoded as submesh #2 now, need to generalize for more other bacteria models)
                        % if col pt is on bacterial tail and no slip, with unknown omega, and doing trick
                        
                        ind = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 7; % coeff for unknown omega for tail rotation relative to body
                        
                        BI_coll(:,ind) = BI_coll(:,ind) + alpha_factor*crossprod(motor_orientation , r_col);
                        %+ alpha_factor*[motor_orientation(2)*r_col(3) - motor_orientation(3)*r_col(2); ...
                        %                                 motor_orientation(3)*r_col(1) - motor_orientation(1)*r_col(3); ...
                        %                                 motor_orientation(1)*r_col(2) - motor_orientation(2)*r_col(1)];
                        % positive sign since the usual negative for this term is negated due to moving it from the RHS into the matrix
                        
                    end
                    
            end
            
        case 2 % free slip, u_c is unknown, move to other side of BIE and add a -1*alpha_factor to its matrix coeff  (alpha_factor is essentially solid angle)
            % note this is same for resistance vs mobility problems
            inds = matrix_props.n_collocation*3 + 3 * (global_u_ind - 1) + [1 2 3]; % corresponding to inds for 3 components of unknown u
            BI_coll(1,inds(1)) = BI_coll(1,inds(1)) - alpha_factor; 
            BI_coll(2,inds(2)) = BI_coll(2,inds(2)) - alpha_factor; 
            BI_coll(3,inds(3)) = BI_coll(3,inds(3)) - alpha_factor;
            
    end
    
    if assembly_input.performance.rigid_body_matrix_rotation && firstrun
        %             RHS_BIE0(global_node_ind,:)           = RHS_coll0(1,:); % RHS_coll = 3 x flowcases , RHS = n_collocation*3 x flowcases
        %             RHS_BIE0(global_node_ind + y_start,:) = RHS_coll0(2,:);
        %             RHS_BIE0(global_node_ind + z_start,:) = RHS_coll0(3,:);
        
        RHS_BIE0_temp(global_node_ind,:,:) = RHS_coll0; %now need x y z as 2rd dim in RHS_BIE0, so no need to transpose anymore
    end
    
    %         RHS_BIE(global_node_ind,:) = RHS_coll(1,:); % RHS_coll = 3 x flowcases , RHS = n_collocation*3 x flowcases
    %         RHS_BIE(global_node_ind + y_start,:) = RHS_coll(2,:);
    %         RHS_BIE(global_node_ind + z_start,:) = RHS_coll(3,:);
    
    RHS_BIE(global_node_ind,:,:) = RHS_coll;
    %             switch assembly_input.BC_type(mesh_ind_coll)
    %                 case 1
    %                     A_temp = BI_coll; % A_a_temp = A_x_temp;  A_b_temp = A_y_temp;  A_c_temp = A_z_temp;
    %                 case 2
    %
    %                     normal = Mesh(mesh_ind_coll).normals(local_vert,:); % unit normal vector at current vert
    %                     tangents = squeeze(Mesh(mesh_ind_coll).tangents(local_vert,:,:)); % 2 orthogonal unit tangent vectors s1, s2 at current vert
    %
    %                     A_temp(1,:) = BI_coll(1,:) * normal(1) + BI_coll(2,:) * normal(2) + BI_coll(3,:) * normal(3);  % u.n at current vert
    %                     % already initialized as zeros  A_b_temp = zeros(size(A_a_temp);  A_c_temp = zeros(size(A_a_temp);
    %                     traction_inds = find( vertind_cols == col_i ); %the 3 consecutive column inds corresponding to col_i vert
    %                     A_temp(2,traction_inds) = tangents(:,1); % 2nd constraint for this vertex is that f.s1 = 0
    %                     A_temp(3,traction_inds) = tangents(:,2);  %3rd constraint for this vertex is that f.s2 = 0
    %             end
    
    
    %%
    switch n_splits
        case 0 % A_BIE is small enough to not need any workaround for stupid Coder array size limitation
            
            %                 A_BIE_1(global_node_ind,:) = BI_coll(1,:); % need to transpose since BI_coll has 3 rows, n_col*3 cols and we want 3 pages, n_col*3 cols in A_BIE_1
            %                 A_BIE_1(global_node_ind + y_start,:) = BI_coll(2,:);
            %                 A_BIE_1(global_node_ind + z_start,:) = BI_coll(3,:);
            
            A_BIE_1(global_node_ind,:,:) = BI_coll;
            
            if assembly_input.performance.rigid_body_matrix_rotation && firstrun
                %                     A_BIE0_1_temp(global_node_ind,:) = BI_coll0(1,:); % need to transpose since BI_coll has 3 rows, n_col*3 cols and we want 3 pages, n_col*3 cols in A_BIE_1
                %                     A_BIE0_1_temp(global_node_ind + y_start,:) = BI_coll0(2,:);
                %                     A_BIE0_1_temp(global_node_ind + z_start,:) = BI_coll0(3,:);
                
                A_BIE0_1_temp(global_node_ind,:,:) = BI_coll0;
            end
        case 1  % split A_BIE into 2 parts column-wise
            %                 n = ceil(matrix_props.n_col*3 / 2);
            N =  3*round( ceil(matrix_props.n_cols / 2) / 3);
            %                 A_BIE_1(global_node_ind,:) = BI_coll(1,1:n);  A_BIE_1(global_node_ind + y_start,:) = BI_coll(2,1:n);  A_BIE_1(global_node_ind + z_start,:) = BI_coll(3,1:n);
            A_BIE_1(global_node_ind,:,:) = BI_coll(:,1:N);
            %                 A_BIE_2(global_node_ind,:) = BI_coll(1,n+1 : end);  A_BIE_2(global_node_ind + y_start,:) = BI_coll(2,n+1 : end);  A_BIE_2(global_node_ind + z_start,:) = BI_coll(3,n+1 : end);
            %              A_BIE_2 = zeros(matrix_props.n_collocation,3,matrix_props.n_cols - size(A_BIE_1,2));
            %              assert(isequal(size(A_BIE_2), [matrix_props.n_collocation,3,matrix_props.n_cols - size(A_BIE_1,2)]));
            A_BIE_2(global_node_ind,:,:) = BI_coll(:,N+1 : end);
            %              A_BIE_2(global_node_ind,:,:) = BI_coll(:,1:n);
            if assembly_input.performance.rigid_body_matrix_rotation && firstrun
                %                     A_BIE0_1(global_node_ind,:) = BI_coll0(1,1:n);  A_BIE0_1(global_node_ind + y_start,:) = BI_coll0(2,1:n);  A_BIE0_1(global_node_ind + z_start,:) = BI_coll0(3,1:n);
                A_BIE0_1_temp(global_node_ind,:,:) = BI_coll0(:,1:N);
                %                     A_BIE0_2(global_node_ind,:) = BI_coll0(1,n+1 : end);  A_BIE0_2(global_node_ind + y_start,:) = BI_coll0(2,n+1 : end);  A_BIE0_2(global_node_ind + z_start,:) = BI_coll0(3,n+1 : end);
                A_BIE0_2_temp(global_node_ind,:,:) = BI_coll0(:,N+1 : end);
            end
        otherwise  %split A_BIE into 4 parts column-wise, or watch code fail if more than 4 splits needed
            %                 n = ceil(matrix_props.n_col*3 / 4);
            N = 3*round(ceil(matrix_props.n_cols / 4)/3);
            %                 A_BIE_1(global_node_ind,:) = BI_coll(1,1:n);  A_BIE_1(global_node_ind + y_start,:) = BI_coll(2,1:n);  A_BIE_1(global_node_ind + z_start,:) = BI_coll(3,1:n);
            A_BIE_1(global_node_ind,:,:) = BI_coll(:,1:N);
            %                 A_BIE_2(global_node_ind,:) = BI_coll(1,n+1 : n * 2);  A_BIE_2(global_node_ind + y_start,:) = BI_coll(2,n+1 : n * 2);  A_BIE_2(global_node_ind + z_start,:) = BI_coll(3,n+1 : n * 2);
            A_BIE_2(global_node_ind,:,:) = BI_coll(:,N+1 : N * 2);
            %                 A_BIE_3(global_node_ind,:) = BI_coll(1,n*2+1 : n * 3);  A_BIE_3(global_node_ind + y_start,:) = BI_coll(2,n*2+1 : n * 3);  A_BIE_3(global_node_ind + z_start,:) = BI_coll(3,n*2+1 : n * 3);
            A_BIE_3(global_node_ind,:,:) = BI_coll(:,N*2+1 : N * 3);
            %                 A_BIE_4(global_node_ind,:) = BI_coll(1,n*3+1 : end);  A_BIE_4(global_node_ind + y_start,:) = BI_coll(2,n*3+1 : end);  A_BIE_4(global_node_ind + z_start,:) = BI_coll(3,n*3+1 : end);
            A_BIE_4(global_node_ind,:,:) = BI_coll(:,N*3+1 : end);
            if assembly_input.performance.rigid_body_matrix_rotation && firstrun
                %                     A_BIE0_1(global_node_ind,:) = BI_coll0(1,1:n);  A_BIE0_1(global_node_ind + y_start,:) = BI_coll0(2,1:n);  A_BIE0_1(global_node_ind + z_start,:) = BI_coll0(3,1:n);
                A_BIE0_1_temp(global_node_ind,:,:) = BI_coll0(:,1:N);
                %                     A_BIE0_2(global_node_ind,:) = BI_coll0(1,n+1 : n * 2);  A_BIE0_2(global_node_ind + y_start,:) = BI_coll0(2,n+1 : n * 2);  A_BIE0_2(global_node_ind + z_start,:) = BI_coll0(3,n+1 : n * 2);
                A_BIE0_2_temp(global_node_ind,:,:) = BI_coll0(:,N+1 : N * 2);
                %                     A_BIE0_3(global_node_ind,:) = BI_coll0(1,n*2+1 : n * 3);  A_BIE0_3(global_node_ind + y_start,:) = BI_coll0(2,n*2+1 : n * 3);  A_BIE0_3(global_node_ind + z_start,:) = BI_coll0(3,n*2+1 : n * 3);
                A_BIE0_3_temp(global_node_ind,:,:) = BI_coll0(:,N*2+1 : N * 3);
                %                     A_BIE0_4(global_node_ind,:) = BI_coll0(1,n*3+1 : end);  A_BIE0_4(global_node_ind + y_start,:) = BI_coll0(2,n*3+1 : end);  A_BIE0_4(global_node_ind + z_start,:) = BI_coll0(3,n*3+1 : end);
                A_BIE0_4_temp(global_node_ind,:,:) = BI_coll0(:,N*3+1 : end);
            end
    end
    
    %%
    
    if assembly_input.performance.debug_mode
%         numels = matrix_props.n_col * matrix_props.n_col*3;
        %                 numels =  matrix_props.n_col * tot_elems;
        
%         if numels < assembly_input.performance.numels_max
            
            %             abs_error_1(col_i, :) = abs_error_temp(1,:);  %the largest abs error over all 54 component integrals
            %             rel_error_1(col_i, :) = rel_error_temp(1,:);  %the largest rel error over all 54 component integrals
            %             fun_evals_1(col_i, :) = fun_evals_temp(1,:);
            %             flag_1(col_i, :) = flag_temp(1,:);
            %
            %             min_value_1(col_i, :) = min_value_temp(1,:);  %smallest (absolute) component integral
            %             max_value_1(col_i, :) = max_value_temp(1,:);  %largest (absolute) component integral
            
            %                    time_1(col_i,:) = time_temp(1,:);  %time spent for each element iter
            
            % outputting all these bastards will be truly painful,
            % so leaving unfinished for now unless it really
            % becomes necessary
%         elseif ceil(numels / 2)  < assembly_input.performance.numels_max
            
%         else
            
%         end
        
    end % debug mode
    
    % can display timing debug info, but only works for non-mexed version
    %         task = getCurrentTask;
    %         ID = task.ID;
    %         disp(['Worker = ',num2str(ID),'      ','col_i = ',num2str(col_i),'   took ',num2str(toc(innertic))]);
    
    
end   %parfor over collocation points %%%%%%%%%%%%%%%%%%%%%%%%%    parfor      %%%%%%%%%%%%%%%%%%%%%%%%%

if assembly_input.performance.rigid_body_matrix_rotation && firstrun
    
    % avoid parfor warning about persistant variable assignments within loop
    %            RHS_BIE0 = reshape( permute(RHS_BIE0_temp,[1 3 2]) , matrix_props.n_collocation*3,matrix_props.n_RHS);
    RHS_BIE0 = reshape( RHS_BIE0_temp , matrix_props.n_collocation*3,matrix_props.n_RHS);
    clear RHS_BIE0_temp;  % check if this performs better than not bothering to clear (prolly not really this, more the A matrices next)
    
    % While reshape doesn't create a copy in memory, permute apparently does, so this isn't the greatest thing to be doing.  But there aren't many
    % alternatives.  If we had separate variables for x, y, z BIE components, we'd still have to combine them into one matrix eventually, requiring
    % a copy in memory.
    % OK, I rejiggered the whole thing to weirdly(?) have x, y, z as the 2nd dim so that we can avoid the permute here and thus any copy of the arrays,
    % at least at this point
    
    % after the reshape, matrix organization is all x, all y, all z BIE components down the rows and usual column arrangement
    
    A_BIE0_1 = reshape( A_BIE0_1_temp , matrix_props.n_collocation*3,size(A_BIE0_1_temp,3));
    clear A_BIE0_1_temp;  % check if this performs better than not bothering
    if ~isempty(A_BIE0_2_temp)
        A_BIE0_2 = reshape( A_BIE0_2_temp , matrix_props.n_collocation*3,size(A_BIE0_2_temp,3));
        clear A_BIE0_2_temp;  % check if this performs better than not bothering
    end
    if ~isempty(A_BIE0_3_temp)
        A_BIE0_3 = reshape( A_BIE0_3_temp , matrix_props.n_collocation*3,size(A_BIE0_3_temp,3));
        clear A_BIE0_3_temp;  % check if this performs better than not bothering
    end
    if ~isempty(A_BIE0_4_temp)
        A_BIE0_4 = reshape( A_BIE0_4_temp , matrix_props.n_collocation*3,size(A_BIE0_4_temp,3));
        clear A_BIE0_4_temp;  % check if this performs better than not bothering
    end
    
end

RHS_BIE = reshape( RHS_BIE , matrix_props.n_collocation*3,matrix_props.n_RHS);

A_BIE_1 = reshape( A_BIE_1 , matrix_props.n_collocation*3,size(A_BIE_1,3));
if ~isempty(A_BIE_2)
    A_BIE_2 = reshape( A_BIE_2 , matrix_props.n_collocation*3,size(A_BIE_2,3));
end
if ~isempty(A_BIE_3)
    A_BIE_3 = reshape( A_BIE_3 , matrix_props.n_collocation*3,size(A_BIE_3,3));
end
if ~isempty(A_BIE_4)
    A_BIE_4 = reshape( A_BIE_4 , matrix_props.n_collocation*3,size(A_BIE_4,3));
end


if assembly_input.performance.rigid_body_matrix_rotation && ~firstrun
    
    
    
    switch assembly_input.problemtype
        case "resistance"
            n_tensors = (matrix_props.n_collocation + matrix_props.n_unknown_u);
        case "mobility"
            n_tensors = (matrix_props.n_collocation + matrix_props.n_unknown_u) + 2; % not sure of correct term for this.  this is the
            % number of 3 x 3 blocks in the matrix that need to be rotated (for each coll pt, we get 3 x 3 coeffs for f, free slip u, U, Omega - each 3 x 3
            % group is a tensor??
    end
    
    %%
    
    DiagM = zeros(3*matrix_props.n_collocation,5); % will contain the 5 diagonals of the final rotation matrix that multiplies the majority of A
    c = 0;
    for d = -2:2 % positions of the 5 diagonals
        c = c + 1;
        diag_entries = diag(rotmat,d); %all entries along diagonal d of the standard rotation matrix
        Diag = NaN(matrix_props.n_collocation * length(diag_entries) , 1); % need to copy each diagonal entry n_tensors times, e.g. # of collocation pts for simple no slip resistance problem
        s = 1;
        for dd = 1:length(diag_entries)
            Diag(s:(s + matrix_props.n_collocation - 1)) = repmat( diag_entries(dd) , matrix_props.n_collocation,1); % fill in the replicated value for each diagonal entry
            s = s + matrix_props.n_collocation;
            
        end
        
        if d <= 0 % this is below the main diagonal, spdiags uses the first values (since entire vector is longer than subdiagonal)
            DiagM(1:matrix_props.n_collocation*length(diag_entries),c) = Diag;
        else % this is above the main diagonal, spdiags uses the last values (since entire vector is longer than superdiagonal)
            DiagM(end - matrix_props.n_collocation*length(diag_entries) + 1:end,c) = Diag;
        end
    end
    
    R = spdiags( DiagM , (-2:2)*matrix_props.n_collocation , matrix_props.n_collocation*3,matrix_props.n_collocation*3);  % sparse with 5 diagonals corresponding to the 5 diagonals of rotmat
    
    R_T = kron(speye(n_tensors),rotmat');  % This is not actually the transpose of R.  Normally, one rotates a tensor M with R A R' but this doesn't work here because of my weird layout of A
    % (all x, all y, all z for rows of BIE but xyz xyz xyz xyz ... across columns for the components of traction, free slip u, U, Omega).  So we get the
    % correct result by defining R' like this.
    %%
    n_non_tensors = (matrix_props.n_cols - n_tensors*3); %equals at most 1 currently, only part of A not to include would be final column for unknown tail rotation omega, if present
    if numel_A_BIE < assembly_input.performance.numels_max  % no splits needed
        
        A_BIE_1(:,1:end - n_non_tensors) = A_BIE_1(:,1:end - n_non_tensors) + R * A_BIE0_1(:,1:end - n_non_tensors) * R_T;
        if assembly_input.rotating_flagellum && strcmp(assembly_input.Tail.motorBC,'torque')
            A_BIE_1(:,n_tensors*3+1) = A_BIE_1(:,n_tensors*3+1) + R * A_BIE0_1(:,n_tensors*3+1); % column for unknown tail rotation omega.  since these are just vector unknowns *(not 3 x 3 tensors), don't need R_T
        end
    elseif ceil(numel_A_BIE / 2)  < assembly_input.performance.numels_max  %split A_BIE into two parts
        A_BIE_1 = A_BIE_1 +  R * A_BIE0_1 * R_T(1:size(A_BIE0_1,2),1:size(A_BIE0_1,2));
        A_BIE_2(:,1:end - n_non_tensors) = A_BIE_2(:,1:end - n_non_tensors) +  R * A_BIE0_2(:,1:end - n_non_tensors) * R_T(1:size(A_BIE0_2,2) - n_non_tensors, 1:size(A_BIE0_2,2) - n_non_tensors);
        if assembly_input.rotating_flagellum && strcmp(assembly_input.Tail.motorBC,'torque')
            ind = n_tensors*3+1 - size(A_BIE0_1,2);
            A_BIE_2(:,ind) = A_BIE_2(:,ind) +  R * A_BIE0_2(:,ind);
        end
    else %ceil(numel_A_BIE / 4)  < assembly_input.performance.numels_max  % split A_BIE into 4 parts
        A_BIE_1 = A_BIE_1 + R * A_BIE0_1 * R_T(1:size(A_BIE0_1,2),1:size(A_BIE0_1,2));
        A_BIE_2 = A_BIE_2 +  R * A_BIE0_2 * R_T(1:size(A_BIE0_2,2),1:size(A_BIE0_2,2));
        A_BIE_3 = A_BIE_3 +  R * A_BIE0_3 * R_T(1:size(A_BIE0_3,2),1:size(A_BIE0_3,2));
        A_BIE_4(:,1:end - n_non_tensors) = A_BIE_4(:,1:end - n_non_tensors) + R * A_BIE0_4(:,1:end - n_non_tensors) * R_T(1:size(A_BIE0_4,2) - n_non_tensors, 1:size(A_BIE0_4,2) - n_non_tensors);
        if assembly_input.rotating_flagellum && strcmp(assembly_input.Tail.motorBC,'torque')
            ind = n_tensors*3+1 - size(A_BIE0_1,2) - size(A_BIE0_2,2) - size(A_BIE0_3,2);
            A_BIE_4(:,ind) = A_BIE_4(:,ind) + R * A_BIE0_4(:,ind);
        end
    end
    
    RHS_BIE = RHS_BIE + R * RHS_BIE0; %this rotates each rotatable x y z block of the original RHS
end



%     if assembly_input.performance.verbose
%         disp(['Matrix assembly for ',Mesh(mesh_ind_coll).name,' collocation pts took ',num2str(toc(mesh_loop_tic))]);
%     end

% end %for over submeshes

if assembly_input.performance.verbose
    disp(['Entire matrix assembly loop took ',num2str(toc)]);
end
tic;

if assembly_input.performance.debug_mode
    %             Relerror = abs(Abserror ./ Integral);
    %             relcontrol = sum(Relerror(:) <= assembly_input.accuracy.integration_tol.traction.reltol & Abserror(:) >= assembly_input.accuracy.integration_tol.traction.abstol) ./ numel(Abserror);
    %             abscontrol = sum(Relerror(:) >= assembly_input.accuracy.integration_tol.traction.reltol & Abserror(:) <= assembly_input.accuracy.integration_tol.traction.abstol) ./ numel(Abserror);
    %             bothcontrol = sum(Relerror(:) <= assembly_input.accuracy.integration_tol.traction.reltol & Abserror(:) <= assembly_input.accuracy.integration_tol.traction.abstol) ./ numel(Abserror);
    
    
    
    %     debug_info.abs_error = abs_error_1;  abs_error_1 = [];
    %     debug_info.rel_error = rel_error_1;  rel_error_1 = [];
    %     debug_info.fun_evals = fun_evals_1;  fun_evals_1 = [];
    %     debug_info.flag      = flag_1;       flag_1 = [];
    %     debug_info.min_value = min_value_1;  min_value_1 = [];
    %     debug_info.max_value = max_value_1;  max_value_1 = [];
    %         debug_info.time = time_1;            time_1 = [];
    
end

% disp(['Debug info and Ax, Ay, Az building took ',num2str(toc)]);



