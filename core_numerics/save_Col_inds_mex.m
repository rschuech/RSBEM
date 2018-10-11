function Col_inds = save_Col_inds_mex(Mesh,nthreads)
%Col_inds is used to map x,y,z traction components at each element vertex to column indices
%of the A matrix

%Col_inds is n_elem x 18 where each row (for each elem) is of the form
% phi1x phi2x phi3x phi4x phi5x phi6x   phi1y phi2y phi3y phi4y phi5y phi6y   phi1z phi2z phi3z phi4z phi5z phi6z

%Col_inds is defined using *global* indices, not local indices

% tot_elems = 0;  tot_verts = 0;
% for i = 1:length(Mesh)  %mex is too dumb to understand sum([Mesh.n_elem])...
% tot_elems = Mesh(i).n_elem + tot_elems;
% tot_verts = Mesh(i).n_vert + tot_verts;
% end

tot_elems = Mesh(end).indices.glob.unq_bounds.elem(2 );  
tot_verts = Mesh(end).indices.glob.unq_bounds.vert(2 );  

temp = repmat(1:tot_verts,3,1);
vertind_cols = temp(:)';  %defined in relation to global indices

Col_inds = NaN(tot_elems,18);

[~, elem_range] = bounds_from_global_inds(Mesh);

parfor (elem_i = 1:tot_elems, nthreads)  %can be a parfor, but slower than serial for small meshes
  % for elem_i = 1:tot_elems
%     imesh = find(elem_i >= elem_range(1,:) & elem_i <= elem_range(2,:)); %which submesh does this global index go with?
%     imesh = imesh(1);  %stop mex from complaining (we know imesh will always be a scalar)
%     local_elem = elem_i - Mesh(imesh).global_indices.elem.start + 1; %local element index going with current global element index
    
      [i_mesh_elem, local_elem] = global2local(elem_i, Mesh, elem_range,'elem');
    
    local_vertinds = Mesh(i_mesh_elem).elems(local_elem,:);  %local vertex indices
    global_vertinds = Mesh(i_mesh_elem).indices.glob.vert(local_vertinds');  %for some reason Coder wants a transpose
    
    %global_vertinds = local_vertinds + Mesh(i_mesh_elem).global_indices.vert.start - 1;  %convert to global vertex indices
    
    col_inds = NaN(3,6);
    for i = 1:6
        col_inds(:,i) = find(vertind_cols == global_vertinds(i));
    end
    
    col_inds = col_inds';
    col_inds = col_inds(:)';
    
    Col_inds(elem_i,:) = col_inds;
    
end
