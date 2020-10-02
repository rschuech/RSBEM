function [VL_unity, A_force] = compute_force_integral(Mesh, matrix_props, index_mapping,input)

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

[VL_unity] = integrate_unity(Mesh, input);


if nargout == 2
    %assemble matrix that can later be multiplied by the solution f to yield total
    %force on body / bodies
    A_force = zeros(3,matrix_props.n_collocation*3,length(Mesh));
    %     [~, elem_range] = bounds_from_global_inds(Mesh);
    
    %     [i_mesh_elems, ~] = global2local(1:tot_elems, Mesh, elem_range, 'elem');    %compute i_mesh_elem and local_elem, indices to current element's submesh and local element index
    
    for i_mesh = 1:length(Mesh) %keep integrals for each submesh separate
        A_force_temp = zeros(3,matrix_props.n_collocation*3);
        VL_unity_submesh = VL_unity{i_mesh};
        % parfor(elem_i = 1:tot_elems, input.performance.nthreads)  %hell, make it parallel anyway - hopefully worth it for huge runs
        %for some unknown reason, mex of the parfor causes a crash when run.
        % it gets all the way through the i_mesh loop but then crashes right
        % at the end
        
        % need to see if above problem still occurs after code revision 9/9/2020
        
        %         for elem_i = 1:tot_elems  %not parallel but who cares - only need to do this once per run, and these are easy operations
        
        parfor elem_i = 1:Mesh(i_mesh).n_elements
            %             i_mesh_elem = i_mesh_elems(elem_i);
            
            
            %             if i_mesh_elem ~= i_mesh %wrong submesh
            %                 continue
            %             end
            element_node_inds = Mesh(i_mesh).elements(elem_i,:); % local inds for nodes of current element  1 x 6
            element_nodes = Mesh(i_mesh).nodes(element_node_inds,:); % coords of nodes of current element   6 x 3
            element_global_node_inds = index_mapping.local_node2global_node{i_mesh}(element_nodes); % global inds for nodes of current element  6 x 1
            inds = ( 3 * (element_global_node_inds - 1) + [1 2 3] )'; % indices of matrix columns corresponding to x,y,z traction components at global nodes of
            % current element   3 x 6  each column is column inds for x y z traction for a node
            
            %             col_inds = matrix_props.Col_inds(elem_i,1:18);
            % sum of forces = 0, entire object
            VL = VL_unity_submesh(elem_i,:);
            %disp('line 46');
            %VL is a 6 element vector that is the integral of hs * phi for this element
            %         A_force(1, col_inds(1:6))   = A_force(1, col_inds(1:6))   + VL;  %x, y, z components have same coeffs    x
            %         A_force(2, col_inds(7:12))  = A_force(2, col_inds(7:12))  + VL;  %x, y, z components have same coeffs    y
            %         A_force(3, col_inds(13:18)) = A_force(3, col_inds(13:18)) + VL;  %x, y, z components have same coeffs    z
            
            
            
            
            temp = zeros(3,matrix_props.n_collocation*3);
            temp(1,inds(1:6)) = VL; % x component?
            temp(2,inds(7:12)) = VL; % y component?
            temp(3,inds(13:18)) = VL; % z component?
            
            A_force_temp = A_force_temp + temp;
            
        end %elems
        
        A_force(:,:,i_mesh) = A_force_temp;
        
    end
    
end

