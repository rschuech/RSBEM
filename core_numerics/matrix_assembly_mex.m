function [A_BIE, RHS_BIE] = matrix_assembly_mex(Mesh,Network, Repulsion, matrix_props,index_mapping,mesh_node_parameters,assembly_input)

%assembles all entries of A_BIE matrix block for both resistance problem and force-free cases

%if debug_mode is true, debug_info can be output

%produces up to 4 submatrices A_BIE_1 - 4 of full A_BIE to avoid mex error
%due to A_BIE being too big.  As a workaround, we insert A_BIE_1 - 4 into full A in
%matrix_assembly_mex_wrapper, in native Matlab, which has no dumb
%limitation on matrix size.
% BI_parameters = assembly_input;
% BI_parameters.Tail = rmfield(BI_parameters.Tail,'motor_torque'); % don't need this for BIEs, only comes in on RHS of additional constraint eqs later
% % and we want to generalize it as a function of time
% BI_parameters = rmfield(BI_parameters,'repulsion');

BI_parameters.performance.eliminate_DL = assembly_input.performance.eliminate_DL;
BI_parameters.performance.DL_singularity_removal = assembly_input.performance.DL_singularity_removal;
BI_parameters.performance.rigid_body_matrix_rotation = assembly_input.performance.rigid_body_matrix_rotation;
BI_parameters.performance.verbose = assembly_input.performance.verbose;
BI_parameters.problemtype = assembly_input.problemtype;
BI_parameters.rotating_flagellum = assembly_input.rotating_flagellum;
BI_parameters.Tail.motorBC = assembly_input.Tail.motorBC;
BI_parameters.Tail.submesh_index = assembly_input.Tail.submesh_index;
BI_parameters.Tail.motor_orientation = assembly_input.Tail.motor_orientation;
BI_parameters.accuracy.mesh = assembly_input.accuracy.mesh;
BI_parameters.accuracy.network = assembly_input.accuracy.network;
BI_parameters.accuracy.triangle_integration = assembly_input.accuracy.triangle_integration;

BI_parameters.constants = assembly_input.constants;


[refpoint, BI_parameters.Tail.motor_orientation] = get_rotational_references(Mesh, assembly_input);




% temp = index_mapping;
% index_mapping = [];
% index_mapping.local_node2global_node = struct2cell(temp.local_node2global_node);
% index_mapping.global_node2local_node = reshape( struct2cell(temp.global_node2local_node) , [] ,1);
% index_mapping.global_node2global_u = temp.global_node2global_u;
% index_mapping.global_u2global_node = temp.global_u2global_node;

% index_mapping = struct('local_node2global_node',{struct2cell(index_mapping.local_node2global_node)},...
%     'global_node2local_node',{reshape( struct2cell(index_mapping.global_node2local_node) , [] ,1)},...
%     'global_node2global_u',index_mapping.global_node2global_u,...
%     'global_u2global_node',index_mapping.global_u2global_node);

% index_mapping.local_node2global_node = struct2cell(index_mapping.local_node2global_node);
% index_mapping.global_node2local_node = struct2cell(index_mapping.global_node2local_node)';
index_mapping2 = index_mapping;
index_mapping2.local_node2global_node = cell2struct(index_mapping.local_node2global_node,'indices',1);
index_mapping2.global_node2local_node = cell2struct(index_mapping.global_node2local_node,'indices',length(index_mapping.global_node2local_node));

Network2 = Network;
Network2.link_members = cell2struct(Network2.link_members,'indices',length(Network2.link_members));


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

% numel_A_BIE = 3 * matrix_props.n_collocation * matrix_props.n_cols; % numel of matrix block corresponding to BIEs
% we get one vector equation (x y z components) by writing the BIE with x0 at a mesh node, and we do this for every node
% need to revise next comment, going back to x y z as 3rd dim due to parfor woes (loop indexing must be the same everywhere within loop)
% A_BIE will be organized with all x components, then all y, then all z (doing xyz xyz xyz ... is much harder to do with parfor)

% y_start = matrix_props.n_collocation; % row index of first y component of BIEs
% z_start = matrix_props.n_collocation*2;  %parfor doesn't allow putting these directly into indexing expressions...

% assembly_input.performance.rigid_body_matrix_rotation = false;


persistent A_BIE0
if assembly_input.performance.rigid_body_matrix_rotation
    if isempty(A_BIE0)
        firstrun = true;
        %         orientations0 = NaN(3,3,length(Mesh));
        %         for i = 1:length(Mesh)
        %             orientations0(:,:,i) = Mesh(i).orientation;
        %         end
    else
        firstrun = false;
    end
else
    
    firstrun = true;
end


% if numel_A_BIE < assembly_input.performance.numels_max  % no splits needed
%     n_splits = 0;
%     A_BIE = zeros(matrix_props.n_collocation,3,matrix_props.n_cols); % x y z components of BIE as 2nd dim to faciliate parfor (which only allows the same
% extremely restrictive (e.g. no : allowed) indexing expression for A_BIE_1 throughout the parfor loop, so no direct way to fill in x, y, z components as different rows)
% weird to have just 3 columns and n_cols pages as 3rd dim but this avoids a presumably costly permute() to eventually make a 2D matrix (permute makes
% a copy in memory, reshape does not)
%     A_BIE_2 = [];    A_BIE_3 = [];     A_BIE_4 = [];

%make sure we can even initialize the full damn thing before going further
A_BIE = zeros(matrix_props.n_collocation, 3, matrix_props.n_cols); % entire BIE matrix block

% NEED TO FIX BELOW FOR NEW MATRIX SETUP (not really cause I didn't end up changing?)
if assembly_input.performance.rigid_body_matrix_rotation && firstrun
    
    
    %             switch assembly_input.performance.rigid_body_tensor_storage
    %                 case 'sparse'
    %                     A_BIE0_temp = sparse(matrix_props.n_collocation,3,matrix_props.n_cols); % BIE0 signifies the values at the initial orientation, to be
    % rotated at later times/phases instead of recomputed from scratch
    % a temp version is used for these because we need x, y, z as 3rd dim to initially fill in entries in the parfor but later need x, y, z
    % as rows to do the tensor rotations and Coder doesn't allow reassignment of persistent variables (e.g. A_BIE0_1 = reshape(A_BIE0_1) )
    % so we'll do initial creation with A_BIE0_temp and then assign A_BIE0 = reshape(A_BIE0_temp) and then use A_BIE0 in later function
    % calls
    %                 case 'full'
    A_BIE0_temp = zeros(matrix_props.n_collocation,3,matrix_props.n_cols); % sparse ND matrices not supported in regular Matlab so will have
    % to store 2 full copies of A temporarily to do this, though we can change it to sparse once we reshape it to 2D and then leave it as
    % sparse for all future timesteps
    %             end
    
    
end



% if matrix_props.n_unknown_u > 0 % some verts have a free slip BC
% u_cols = zeros(matrix_props.n_col, matrix_props.n_unknown_u*3 ,3); % extra columns for coeffs corresponding to unknown u components in DL integrals
% u_rows = zeros(matrix_props.n_unknown_u*3 , matrix_props.n_cols);
% this is now part of BI_coll variable, lumped with traction

% end

RHS_BIE = zeros( matrix_props.n_collocation, 3, matrix_props.n_RHS); % add extra fake rows to match setup of A_BIE
% need to revise for x y z organization% to match A_BIE, organized as all x, all y, all z, and 2nd dim is for possible multiple RHS

persistent RHS_BIE0
if assembly_input.performance.rigid_body_matrix_rotation
    if firstrun
        RHS_BIE0_temp = zeros(matrix_props.n_collocation,3,matrix_props.n_RHS);
    end
else
    RHS_BIE0_temp = [];
end


if assembly_input.performance.verbose
    disp('initialized A_BIE');
end



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





% BI_parameters.constants.refpoint = refpoint;
% BI_parameters.x0_location = "on_mesh";
% BI_parameters.performance.rigid_body_matrix_rotation = assembly_input.performance.rigid_body_matrix_rotation;
% BI_parameters.accuracy.mesh.ignore_interaction = assembly_input.accuracy.mesh.ignore_interaction;
% BI_parameters.accuracy.mesh.eps2 = assembly_input.accuracy.mesh.eps2;
% BI_parameters.accuracy.network.eps2

x0_location = "on_mesh";
BI_parameters.constants.refpoint = refpoint;

% coder.varsize('reference_nodes',[2 3 Inf],[true true true]); %for mex rules
% rng(0);
% rand_coll_inds = randperm(matrix_props.n_collocation);
 rand_coll_inds = 1:matrix_props.n_collocation;

parfor (rand_global_node_ind = 1:matrix_props.n_collocation , assembly_input.performance.nthreads)  %rows of A_BIE in sets of 3, i.e. x y z components of BIE for each node
   
    global_node_ind = rand_coll_inds(rand_global_node_ind);
    
    % this loops over collocation points
    %         for global_node_ind = 1:matrix_props.n_collocation    %randperm(matrix_props.n_collocation)
    % instead of randomizing indices of actual meshes, why not introduce a vector that is a random permutation of 1:n_collocation and then here, use
    % randvec(global_node_ind) to get the actual collocation pt index that we get vert coords etc for.  This randomizes order of loop over collocation pts but
    % keeps meshes in original (presumably somewhat logical) node order.  After matrix solve, can un-randomize order of solution vector so that
    % traction, free slip velocity, etc match the original meshes.
    
    % need to use below again in separate code that calcs u at network nodes
    %     switch BI_parameters.x0_location
    %         case "on_mesh"
    %             temp = unique( index_mapping.global_node2local_node{global_node_ind}(:,[1 2]), 'rows' ); % submeshes and local node inds for this global node
    %             coll_local = struct('submeshes',temp(:,1),'nodes',temp(:,2));
    %         case "off_boundary"
    %             coll_local = struct('submeshes',NaN(0,1),'nodes',NaN(0,2));
    %     end
    
    temp = unique( index_mapping.global_node2local_node{global_node_ind}(:,[1 2]), 'rows' ); % submeshes and local node inds for this global node
    x0_parameters = struct('x0',Mesh(temp(1,1)).nodes(temp(1,2),:)',  'submeshes',temp(:,1),      'nodes',temp(:,2));
    x0_parameters.location = x0_location;
    
    
    
    x0_parameters.r = Mesh(temp(1,1)).nodes(temp(1,2),:)' - BI_parameters.constants.refpoint; % vector from refpoint to collocation pt
    
    
    switch x0_location
        case "on_mesh"
            x0_parameters.BC_type = mesh_node_parameters.BC_type(global_node_ind);
            % index for the unknown velocity at the current collocation point,
            % between 1 - n where n is number of free slip nodes
            x0_parameters.global_u_ind = index_mapping.global_node2global_u( global_node_ind );
            x0_parameters.u = reshape(  squeeze(Mesh(temp(1,1)).u(temp(1,2),:,:)), 3 , []  ); % needed if including DL  3 x flowcases
            
        case "off_mesh"
            x0_parameters.BC_type = NaN;
            x0_parameters.global_u_ind = NaN;
            x0_parameters.u = NaN(3,1); % u_coll will never be used
    end
    
    
    % reshape needed in case there's only one flowcase, to make sure
    % arbitrarily look up u_coll for first submesh/local node this global node appears in - u_coll should be the same for all appearances
    
    
    
    [BI_coll, RHS_coll, BI_coll0, RHS_coll0] = boundary_integrals_mexed(x0_parameters, global_node_ind, Mesh, mesh_node_parameters, Network2, Repulsion, matrix_props, index_mapping2, BI_parameters ,firstrun);
    
    
    [full_DL_submeshes_containing_collocation_pt, inds] = setdiff( x0_parameters.submeshes , find(assembly_input.performance.eliminate_DL | assembly_input.performance.DL_singularity_removal == 1) ); % submeshes that this collocation pt is on for which we are not eliminating the DL or using trick
    % note, eliminate_DL should always be either true or false while DL_singularity_removal can be NaN for bodies that should always have DL eliminated
    % e.g. sheets
    alphas = NaN(length(full_DL_submeshes_containing_collocation_pt),1);
    for n = 1:length(full_DL_submeshes_containing_collocation_pt)
        alphas(n) = Mesh( full_DL_submeshes_containing_collocation_pt(n) ).solid_angle(x0_parameters.nodes(inds(n)));
    end
    % note that we may still be doing full (or singularity removal trick) DL calcs on submeshes NOT containing collocation pt
    alpha_factor = sum(alphas) + (1 - length(full_DL_submeshes_containing_collocation_pt))*4*pi; % this is what multiplies the collocation velocity u(x_0) in the BIE; in the typical case
    % of u(x_0) being a member of just one body, alpha_factor = alpha.  For e.g. u(x_0) on dino body and transverse, since transverse alpha = 4*pi, we
    % still get alpha of just the body here.
    
    %     if assembly_input.performance.eliminate_DL(coll_local.submesh) || assembly_input.performance.DL_singularity_removal
    %         alpha = 4*pi; % whether this mesh is a sheet, a rigid body, or a deforming body (that we don't need correct traction for), eliminating the
    %         %DL results in solid angle effectively becoming 4*pi regardless of the actual geometric value
    %         %keeping the DL integral but using the singularity removal identity for the mesh containing the col pt also results in alpha effectively becoming 4*pi
    %     else
    %         alpha = Mesh(coll_local.submesh).solid_angle(coll_local.node);  % solid angle of surface at current collocation pt:  2*pi for smooth, closer to zero for a valley, closer to 4*pi for a ridge or peak
    %     end
    
    
    % currently redoing below contributions even for rigid body coll pts + elements since I'm not sure if/how the rotations would work for these.
    % the cost of redoing them is probably very small though
    
    % add matrix and/or RHS contributions from the alpha_factor*u_c = side of the BIEs, whether alpha_factor is typically the actual solid angle at the coll pt or set to
    % 4 pi due to eliminating the DL or using the singularity removal identity on the mesh containing the coll pt
    % note that if using singularity removal identity on any meshes NOT containing the coll pt, then the last term in the identity eq is zero for each of these DL integrals
    % (it is (alpha - 4*pi)u_c and that alpha is 4 pi since the coll point is not in/on that mesh - note, not same alpha_factor as the LHS of the BIE) so we
    % don't need to do anything about it here
    
    switch x0_parameters.BC_type
        case 1 % no slip
            % RHS_elem is 3 x flowcases
            RHS_coll = RHS_coll + alpha_factor*x0_parameters.u; % contribution of known u_r    u_col is 3 x flowcases
            
            switch assembly_input.problemtype
                
                %                         case 'resistance' % known u_c goes on RHS
                
                case "mobility" % u_c = U + (Omega x r) + u_r
                    inds = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + (1:6);
                    % subtract to move unknowns into matrix side of BIE
                    BI_coll(1,inds ) = BI_coll(1,inds ) - alpha_factor*[1 0 0   0          x0_parameters.r(3)   -x0_parameters.r(2)]; % coeffs from U and (Omega x r)
                    BI_coll(2,inds ) = BI_coll(2,inds ) - alpha_factor*[0 1 0  -x0_parameters.r(3)   0           x0_parameters.r(1)]; % coeffs from U and (Omega x r)
                    BI_coll(3,inds ) = BI_coll(3,inds ) - alpha_factor*[0 0 1   x0_parameters.r(2)  -x0_parameters.r(1)    0]; % coeffs from U and (Omega x r)
                    
                    if assembly_input.rotating_flagellum && assembly_input.Tail.motorBC == "torque" && ismember( assembly_input.Tail.submesh_index , x0_parameters.submeshes)
                        % (tail hardcoded as submesh #2 now, need to generalize for more other bacteria models)
                        % if col pt is on bacterial tail and no slip, with unknown omega, and doing trick
                        
                        ind = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 7; % coeff for unknown omega for tail rotation relative to body
                        
                        BI_coll(:,ind) = BI_coll(:,ind) - alpha_factor*crossprod(BI_parameters.Tail.motor_orientation , x0_parameters.r);
                        %+ alpha_factor*[motor_orientation(2)*r_col(3) - motor_orientation(3)*r_col(2); ...
                        %                                 motor_orientation(3)*r_col(1) - motor_orientation(1)*r_col(3); ...
                        %                                 motor_orientation(1)*r_col(2) - motor_orientation(2)*r_col(1)];
                        % WRONG I THINK:  positive sign since the usual negative for this term is negated due to moving it from the RHS into the matrix
                        % I think there's a mistake in the PNAS SI, this term should normally have a postive sign, but it gets negated when moving it from
                        % u-side of BIE into the matrix side.  Having a negative here matches old code and appears to work as far as particle tracking
                        % tests.
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
        
        RHS_BIE0_temp(rand_global_node_ind,:,:) = RHS_coll0; %now need x y z as 2rd dim in RHS_BIE0, so no need to transpose anymore
    end
    
    
    RHS_BIE(rand_global_node_ind,:,:) = RHS_coll;
    
    
    
    
    A_BIE(rand_global_node_ind,:,:) = BI_coll;
    
    if assembly_input.performance.rigid_body_matrix_rotation && firstrun
        
        A_BIE0_temp(rand_global_node_ind,:,:) = BI_coll0;
    end
    
    
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
    switch assembly_input.performance.rigid_body_tensor_storage
        case 'sparse'
            A_BIE0 = sparse(reshape( A_BIE0_temp , matrix_props.n_collocation*3,size(A_BIE0_temp,3)));
        case 'full'
            A_BIE0 = reshape( A_BIE0_temp , matrix_props.n_collocation*3,size(A_BIE0_temp,3));
    end
    
    clear A_BIE0_temp;  % check if this performs better than not bothering
    
end


RHS_BIE = reshape( RHS_BIE , matrix_props.n_collocation*3, matrix_props.n_RHS);
% RHS_BIE(matrix_props.n_rows + 1 : end,:) = [];  % delete extra fake rows, if there are any

A_BIE = reshape( A_BIE , matrix_props.n_collocation*3, size(A_BIE,3)); % has extra fake rows beyond the constraint eqs to be added
% A_BIE(matrix_props.n_rows + 1 : end, :) = []; % delete extra fake rows, if there are any

% if ~isempty(A_BIE_2)
%     A_BIE_2 = reshape( A_BIE_2 , matrix_props.n_collocation*3,size(A_BIE_2,3));
% end
% if ~isempty(A_BIE_3)
%     A_BIE_3 = reshape( A_BIE_3 , matrix_props.n_collocation*3,size(A_BIE_3,3));
% end
% if ~isempty(A_BIE_4)
%     A_BIE_4 = reshape( A_BIE_4 , matrix_props.n_collocation*3,size(A_BIE_4,3));
% end


if assembly_input.performance.rigid_body_matrix_rotation && ~firstrun
    % rotation_matrices = 3 x 3 x length(Mesh) describe how to rotate initial meshes to get to current orientation
    %     if size(rotation_matrices,3) ~= length(Mesh)
    %         if size(rotation_matrices,3) == 1
    %             rotation_matrices = repmat(rotation_matrices,1,1,length(Mesh)); % expand single rotation matrix, assuming it applies to all submeshes
    %         else
    %             error('Check number of rotation matrices input, should be either 1 or # submeshes');
    %         end
    %     end
    
    first_submesh = @(x) x(1,1);  % x will be a cell in index_mapping.global_node2local_node.  We will arbitrarily take the first submesh listed.
    % If a global node is shared by a few submeshes, each of those submeshes better have the same rotmat or else there are bigger problems....
    %     global_node2submesh = NaN(matrix_props.n_collocation,1);
    %     for i = 1:length(Mesh)
    global_node2submesh = cellfun(first_submesh, index_mapping.global_node2local_node);  % n_collocation x 1 of which rotmat each global node goes with
    % not accounting for randomization of rows of A yet but will below
    
    switch assembly_input.problemtype
        case "resistance"
            n_vector_unknowns = (matrix_props.n_collocation + matrix_props.n_unknown_u);
        case "mobility"
            n_vector_unknowns = (matrix_props.n_collocation + matrix_props.n_unknown_u) + 2; % not sure of correct term for this.  this is the
            % number of 3 x 3 blocks in the matrix that need to be rotated (for each coll pt, we get 3 x 3 coeffs for f, free slip u, U, Omega - each 3 x 3
            % group is a tensor??
    end
    n_scalar_unknowns = (matrix_props.n_cols - n_vector_unknowns*3); %equals at most 1 currently, only part of A not to include would be final column for unknown tail rotation omega, if present
    
    %%
    index_map = [3 NaN NaN; 2 6 NaN; 1 5 9; 4 8 NaN; 7 NaN NaN]; % row 1:5 goes with diag -2:2), col 1:3 goes with diag entry, values are linear indices of rotmat
    diag_length = [1 2 3 2 1];
    DiagM = zeros(3*matrix_props.n_collocation,5); % will contain the 5 diagonals of the final rotation matrix that multiplies the majority of A
    % DiagM is as big (# rows) as the longest, middle diagonal, while the other diagonals are shorter, so we will not fill in this entire matrix
    c = 0;
    for d = -2:2 % positions of the 5 diagonals
        c = c + 1;
        %         diag_entries = diag(rotmat,d); %all entries along diagonal d of the standard rotation matrix
        Diag = NaN(matrix_props.n_collocation * diag_length(c) , 1); % need to copy each rotmat diagonal entry n_collocation times to go with rows of A
        
        s = 1;
        for dd = 1:diag_length(c) % either 1, 2, 3, 2, 1 depending on which diagonal of rotmat we're on
            diag_entries = NaN(matrix_props.n_collocation,1);
            for i = 1:length(Mesh)
                %                 diag_entries = diag( rotation_matrices(:,:,i), d );
                [row,col] = ind2sub([3 3],index_map(c,dd)); % could also just have directly stored row and col inds in index_map
                diag_entries(global_node2submesh(rand_coll_inds) == i) = Mesh(i).rotation_matrix(row,col);
                % here is where we account for randomized rows of A
            end
            %             Diag(s:(s + matrix_props.n_collocation - 1)) = repmat( diag_entries(dd) , matrix_props.n_collocation,1); % fill in the replicated value for each diagonal entry
            Diag(s:(s + matrix_props.n_collocation - 1)) = diag_entries;
            s = s + matrix_props.n_collocation;
            
        end
        
        if d <= 0 % this is below the main diagonal, spdiags uses the first values (since entire vector is longer than subdiagonal)
            DiagM(1:matrix_props.n_collocation*diag_length(c),c) = Diag;
        else % this is above the main diagonal, spdiags uses the last values (since entire vector is longer than superdiagonal)
            DiagM(end - matrix_props.n_collocation*diag_length(c) + 1:end,c) = Diag;
        end
    end
    
    R = spdiags( DiagM , (-2:2)*matrix_props.n_collocation , matrix_props.n_collocation*3,matrix_props.n_collocation*3);  % sparse with 5 diagonals corresponding to the 5 diagonals of rotmat
    
    
    % R_T doesn't need any modification due to randomized rows of A because
    % we take each row of A * each col of R_T.  The order of the columns of
    % A is important here, not the rows, and we haven't changed the order
    % of the columns.
    
    %     R_T = kron(speye(n_vector_unknowns),rotmat');  % This is not actually the transpose of R.  Normally, one rotates a tensor M with R A R' but this doesn't work here because of my weird layout of A
    % (all x, all y, all z for rows of BIE but xyz xyz xyz xyz ... across columns for the components of traction, free slip u, U, Omega).  So we get the
    % correct result by defining R' like this.
    %     R_T = sparse(n_vector_unknowns*3, n_vector_unknowns*3);
    blocks = cell(n_vector_unknowns,1);  blocks(:) = {NaN(3,3)};
    
    inds_collocation = 1:matrix_props.n_collocation;
    inds_unknown_u = (matrix_props.n_collocation + 1) : (matrix_props.n_collocation + matrix_props.n_unknown_u);
    for i = 1:length(Mesh)
        blocks( inds_collocation(global_node2submesh == i) ) = {sparse( Mesh(i).rotation_matrix(:,:)' )};
        blocks( inds_unknown_u( global_node2submesh( index_mapping.global_u2global_node ) == i ) ) = {sparse( Mesh(i).rotation_matrix(:,:)' )};
    end
    
    
    % assume that any blocks to be rotated after the traction and unknown u blocks (i.e. U, Omega) go with rotation matrix for Body
    blocks(matrix_props.n_collocation + matrix_props.n_unknown_u + 1 : end) = {sparse( Mesh([Mesh.name] == "Body").rotation_matrix(:,:)' )};
    
    
    R_T = blkdiag(blocks{:});
    
    
    %     for i = 1:matrix_props.n_collocation
    %         R_T(3*(i-1)+1:3*(i-1)+3,3*(i-1)+1:3*(i-1)+3) = rotation_matrices(:,:,global_node2submesh(i))
    
    %%
    %     if numel_A_BIE < assembly_input.performance.numels_max  % no splits needed
    
    

    
    A_BIE(:,1:end - n_scalar_unknowns) = A_BIE(:,1:end - n_scalar_unknowns) + R * A_BIE0(:,1:end - n_scalar_unknowns) * R_T;
    if assembly_input.rotating_flagellum && strcmp(assembly_input.Tail.motorBC,'torque')
        A_BIE(:,n_vector_unknowns*3+1) = A_BIE(:,n_vector_unknowns*3+1) + R * A_BIE0(:,n_vector_unknowns*3+1); % column for unknown tail rotation omega.  since these are just vector unknowns *(not 3 x 3 tensors), don't need R_T
    end
    
    RHS_BIE(:) = RHS_BIE(:) + R * RHS_BIE0(:); %this rotates each rotatable x y z block of the original RHS
end  % rigid body recycling shortcut and not first run


% firstrun
if assembly_input.performance.verbose
    disp(['Entire matrix assembly loop took ',num2str(toc)]);
end


