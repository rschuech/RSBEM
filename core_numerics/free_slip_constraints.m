function [A_free_slip, RHS_free_slip] = free_slip_constraints(Mesh, matrix_props,index_mapping,node_parameters,assembly_input)



% no_penetration = spalloc(matrix_props.n_unknown_u, matrix_props.n_cols, 10*matrix_props.n_unknown_u); % matrix rows for u.n = 0 (when stated in moving frame)
% no_shear_stress1 = spalloc(matrix_props.n_unknown_u, matrix_props.n_cols, 3*matrix_props.n_unknown_u); % matrix rows for f.s1 = 0 where s1 is a tangent vector
% no_shear_stress2 = spalloc(matrix_props.n_unknown_u, matrix_props.n_cols, 3*matrix_props.n_unknown_u); % matrix rows for f.s2 = 0 where s1 is a tangent vector

RHS_free_slip = zeros(matrix_props.n_unknown_u*3, matrix_props.n_RHS);


A_free_slip = spalloc(matrix_props.n_unknown_u*3, matrix_props.n_cols, 16*matrix_props.n_unknown_u);
% 3 constraint eqs per free slip node:  u.n = 0 (as stated in body frame, no-penetration), f.s1 = 0, f.s2 = 0 where s1, s2 are a tangent vectors
% should be at most 10 non-zero coeffs for u.n = 0 (3 for unknown u, 3 for U, 3 for Omega, 1 for omega) and 3 each for f.s1 = 0, f.s2 = 0, so 16 total per
% free slip node


% for N nodes with unknown free slip u, have N rows of u.n constraint, then N rows of f.s1 = 0, then N rows of f.s2 = 0

% global_u_mesh_inds = NaN(matrix_props.n_unknown_u, 1);
% local_u_inds = NaN(matrix_props.n_unknown_u, 1);
% for global_u_ind = 1:matrix_props.n_unknown_u

%     for i_mesh = 1:length(Mesh)
%         if ismember(global_u_ind,Mesh(i_mesh).indices.glob.unknown_u)


% global_u_inds = [];  global_u_mesh_inds = [];  local_u_inds = [];
% for i_mesh = 1:length(Mesh)
%     temp = Mesh(i_mesh).indices.glob.unknown_u(Mesh(i_mesh).indices.glob.unknown_u ~= 0);
%     global_u_inds = [global_u_inds; temp ];
%     local_u_inds = [local_u_inds; (1:Mesh(i_mesh).n_nodes)'];
%     global_u_mesh_inds = [global_u_mesh_inds; repmat(i_mesh,length(temp),1)];
% end

% [global_u_inds, inds] = unique(global_u_inds);
% global_u_mesh_inds = global_u_mesh_inds(inds);
% local_u_inds = local_u_inds(inds);



[refpoint, motor_orientation] = get_rotational_references(Mesh, assembly_input);

n_global_u = max(index_mapping.global_node2global_u); % number of nodes with unknown u across all submeshes

%       parfor (global_u_ind = global_u_inds , assembly_input.performance.nthreads)  %rows of A_BIE in sets of 3, i.e. x y z components of BIE for each node
for global_u_ind = 1:n_global_u  % parfor would require some temp variables and prolly not worth bothering?
    global_node_ind = find( index_mapping.global_node2global_u == global_u_ind ); % since possibly not all nodes are free-slip, global_u_ind ~= global_node_ind
    submesh_ind = index_mapping.global_node2local_node{global_node_ind}(1,1); % arbitrarily choose the first local instance of this node, since all instances
    % had better have the same u
    local_node_ind = index_mapping.global_node2local_node{global_node_ind}(1,2);
    normal = node_parameters.normals_avg(global_node_ind,:); % unit normal vector
    tangents = squeeze(node_parameters.tangents_avg(global_node_ind,:,:)); % 2 orthogonal unit tangent vectors s1, s2
    inds = matrix_props.n_collocation*3 + 3 * (global_u_ind - 1) + [1 2 3]; % corresponding to inds for 3 components of unknown u
    A_free_slip(global_u_ind,inds) = normal; % u . n  where u is unknown fluid velocity
    RHS_free_slip(global_u_ind,:) = normal * reshape( squeeze(Mesh(submesh_ind).u( local_node_ind , : ,:)) , 3 , []); % n . u = u . n where u is known boundary velocity
    
    switch assembly_input.problemtype
        
        
        case "mobility"
            
            inds = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + (1:6);
            r_col = Mesh(submesh_ind).nodes(local_node_ind,:)' - refpoint; % vector from refpoint to collocation pt
            
            A_free_slip(global_u_ind,inds) = - [ normal,  crossprod(r_col , normal)' ]; % -U . n and -(Omega x r) . n
            % note that Omega x r . n = (r x n) . Omega, so r x n are the coeffs for unknown Omega
            % negative in front of both since we're moving these to the LHS and into the matrix
            
            if assembly_input.rotating_flagellum && strcmp(assembly_input.Tail.motorBC,'torque') && Mesh(submesh_ind).name == "Tail"
                % (assuming one tail for now, need to generalize for multiflagellated)
                % here I'm assuming that the rotating Tail is one completely separate (non coincident) submesh from the others, so that nodes on the Tail are
                % only members of the Tail submesh.  If we have multiple rotating tails or are separating the tail into multiple submeshes, will need to
                % revise this.
                
                ind = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 7; %matrix column for coeff for unknown omega for tail rotation relative to body
                
                A_free_slip(global_u_ind,ind) = normal * crossprod(motor_orientation, r_col)'; % n . (e x r) where e is motor orientation
                % would be negative in front due to moving to LHS but this term already has a negative sign as per definition, see Schuech et al PNAS SI
                
            end
    end
    
    traction_inds = matrix_props.n_collocation*3 + (global_u_ind - 1)*3 + (1:3); %the 3 consecutive matrix column inds corresponding to this node
    
    %     A_constraints([1 2] .* (matrix_props.n_unknown_u + global_u_ind) , traction_inds) = tangents';
    A_free_slip(global_u_ind + matrix_props.n_unknown_u, traction_inds) = tangents(:,1);
    A_free_slip(global_u_ind + matrix_props.n_unknown_u*2, traction_inds) = tangents(:,2);
    % no need to do anything special for RHS here since f.s1 = f.s2 = 0
end




