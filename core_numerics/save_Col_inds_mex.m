function Col_inds = save_Col_inds_mex(Mesh,nthreads)
%Col_inds is used to map x,y,z traction components at each element vertex to column indices
%of the A matrix

%Col_inds is n_elem x 18 where each row (for each elem) is of the form
% phi1x phi2x phi3x phi4x phi5x phi6x   phi1y phi2y phi3y phi4y phi5y phi6y   phi1z phi2z phi3z phi4z phi5z phi6z

% note that this is NOT how the unknowns are ordered in the A matrix; in
% A, they are ordered 1x 1y 1z 2x 2y 2z 3x 3y 3z 4x 4y 4z ... Nx Ny Nz
% in order of the global vert list, all the way to the final global vert
% Col_inds gives the global inds corresponding to the columns of A for each
% elem as we loop over elements and integrate during matrix assembly


%Col_inds is defined using *global* indices, not local indices

tot_elems = Mesh(end).indices.glob.unq_bounds.elem(2 );
tot_verts = Mesh(end).indices.glob.unq_bounds.vert(2 );

temp = repmat(1:tot_verts,3,1); % [1:N; 1:N; 1:N] for x y z components for all verts
vertind_cols = temp(:)';  %defined in relation to global indices  [x y z x y z x y z ...] or [1 1 1 2 2 2 3 3 3 ... N N N] where N is final global vert

Col_inds = NaN(tot_elems,18);

[~, elem_range] = bounds_from_global_inds(Mesh);

parfor (elem_i = 1:tot_elems, nthreads)  %can be a parfor, but slower than serial for small meshes
    % for elem_i = 1:tot_elems
    
    [i_mesh_elem, local_elem] = global2local(elem_i, Mesh, elem_range,'elem'); % get which submesh and which local elem index goes with current global elem
    
    local_vertinds = Mesh(i_mesh_elem).elems(local_elem,:);  %local vertex indices for current element
    global_vertinds = Mesh(i_mesh_elem).indices.glob.vert(local_vertinds');  %for some reason Coder wants a transpose
    
    col_inds = NaN(3,6); % x y z  by  6 verts of current elem
    for i = 1:6 % each vert of current elem
        col_inds(:,i) = find(vertind_cols == global_vertinds(i)); % find inds from vertind_cols (corresponding to x, y, z global inds) matching current vert of current elem 
        % same inds are copied for each of 3 rows for x, y, z components
        % for f at current global vert
    end
    
    col_inds = col_inds'; % now 6 x 3
    col_inds = col_inds(:)'; % [vert_1_x vert_2_x ... vert_6_x   vert_1_y vert_2_y ... vert_6_y   vert_1_z vert_2_z ... vert_6_z]
    
    Col_inds(elem_i,:) = col_inds;  % as described at beginning
    
end
