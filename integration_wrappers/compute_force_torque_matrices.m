function [A_force, A_torque] = compute_force_torque_matrices(Mesh, matrix_props, index_mapping,input)
%%
% tic
%unlike other integration functions, this and compute_torque_integral are
%designed to only take entire Mesh as input, since it must go with solution
%f to be useful and f is presumably for entire Mesh

%A_force is 3 x n_col*3 x length(Mesh)
%for total force on all submeshes, use sum(A_force,3)

% if length(Mesh) > 1 %multiple submeshes input, use global indices
%     tot_elems = Mesh(end).indices.glob.unq_bounds.elem(2 );
% else %only one submesh input, use local indices
%     tot_elems = Mesh(1).n_elem;
% end

[VL_unity] = integration_wrapper(Mesh, 'unity', ...
    input.accuracy.mesh.integration_tol.force.abstol,...
    input.accuracy.mesh.integration_tol.force.reltol,...
    input.accuracy.mesh.integration_tol.force.maxevals,...
    input.performance.nthreads, input.accuracy.triangle_integration); %compute integrals of hS * phi

[refpoint] = get_rotational_references(Mesh, input);

% [VL_unity] = integrate_unity(Mesh, input);
[VL_moments] = integration_wrapper(Mesh, 'moments', ...
    input.accuracy.mesh.integration_tol.torque.abstol,...
    input.accuracy.mesh.integration_tol.torque.reltol,...
    input.accuracy.mesh.integration_tol.torque.maxevals,...
    input.performance.nthreads, input.accuracy.triangle_integration, ...
    refpoint); %compute integrals of hS * phi

% if nargout == 2
    %assemble matrix that can later be multiplied by the solution f to yield total
    %force on body / bodies
    A_force = zeros(3,matrix_props.n_collocation*3,length(Mesh));
    A_torque = zeros(3,matrix_props.n_collocation*3,length(Mesh));
    %     [~, elem_range] = bounds_from_global_inds(Mesh);
    
    %     [i_mesh_elems, ~] = global2local(1:tot_elems, Mesh, elem_range, 'elem');    %compute i_mesh_elem and local_elem, indices to current element's submesh and local element index
    
    for i_mesh = 1:length(Mesh) %keep integrals for each submesh separate
        A_force_temp = zeros(3,matrix_props.n_collocation*3);
        VL_unity_submesh = VL_unity{i_mesh};
        
        A_torque_temp = zeros(3,matrix_props.n_collocation*3);
        VL_moments_submesh = VL_moments{i_mesh};
        % parfor(elem_i = 1:tot_elems, input.performance.nthreads)  %hell, make it parallel anyway - hopefully worth it for huge runs
        %for some unknown reason, mex of the parfor causes a crash when run.
        % it gets all the way through the i_mesh loop but then crashes right
        % at the end
        if any(isnan(VL_moments_submesh(:)))
            stpa
        end
        % need to see if above problem still occurs after code revision 9/9/2020
        
        %         for elem_i = 1:tot_elems  %not parallel but who cares - only need to do this once per run, and these are easy operations
     
        for elem_i = 1:Mesh(i_mesh).n_elements  % this was a parfor but that results in nondeterministic output since
            % it contains a running sum and if the order is random, the
            % roundoff errors will be random
            % seems like there's little difference in run time, if anything
            % for may be faster than parfor for simple spherical swimmer

            element_node_inds = Mesh(i_mesh).elements(elem_i,:); % local inds for nodes of current element  1 x 6
        
            element_global_node_inds = index_mapping.local_node2global_node{i_mesh}(element_node_inds); % global inds for nodes of current element  6 x 1
            inds = ( 3 * (element_global_node_inds - 1) + [1 2 3] ); % indices of matrix columns corresponding to x,y,z traction components at global nodes of
            % current element    6 x 3 each column are the matrix column inds of the 6 nodes, and we have 3 columns for x, y, z components of traction
            
         
            %VL is a 6 element vector that is the integral of hs * phi for this element
            %         A_force(1, col_inds(1:6))   = A_force(1, col_inds(1:6))   + VL;  %x, y, z components have same coeffs    x
            %         A_force(2, col_inds(7:12))  = A_force(2, col_inds(7:12))  + VL;  %x, y, z components have same coeffs    y
            %         A_force(3, col_inds(13:18)) = A_force(3, col_inds(13:18)) + VL;  %x, y, z components have same coeffs    z
            
            
            temp = zeros(3,matrix_props.n_collocation*3);
            % below uses linear indexing, as if we did inds(:) and indexed into that vector
            % the vector would be [6 column inds for x components; 6 column inds for y components; 6 column inds for z components]
            temp(1,inds(1:6)) = VL_unity_submesh(elem_i,:); % x component?
            temp(2,inds(7:12)) = VL_unity_submesh(elem_i,:); % y component?
            temp(3,inds(13:18)) = VL_unity_submesh(elem_i,:); % z component?
            
            A_force_temp = A_force_temp + temp;
            
            temp = zeros(3, matrix_props.n_collocation*3);
            temp(1,inds(13:18)) = VL_moments_submesh(elem_i,7:12);
            temp(1,inds(7:12)) = -VL_moments_submesh(elem_i,13:18);
            temp(2,inds(1:6)) = VL_moments_submesh(elem_i,13:18);
            temp(2,inds(13:18)) = -VL_moments_submesh(elem_i,1:6);
            temp(3,inds(7:12)) = VL_moments_submesh(elem_i,1:6);
            temp(3,inds(1:6)) = -VL_moments_submesh(elem_i,7:12);
            
            A_torque_temp = A_torque_temp + temp;
            % perhaps this can be changed back to parfor with deterministic
            % A_torque_temp and A_force_temp if instead of a cumulative sum
            % we store each value in a big matrix and then add it all up
            % outside the loop
            
        end %elems
        
        A_force(:,:,i_mesh) = A_force_temp;
        A_torque(:,:,i_mesh) = A_torque_temp;
        
    end
    
% end
% toc
