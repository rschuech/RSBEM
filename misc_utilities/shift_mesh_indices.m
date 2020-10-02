function [Mesh, original_indices] = shift_mesh_indices(names,Mesh, original_indices)

ind = find(strcmp('Tail',names));  %Mesh ind going with tail, if we loaded tail
if ~isempty(ind)
    other_inds = setdiff(1:length(names), ind);  % inds for all other submeshes
    %tail indices might overlap body/transverse/wingtip since they're not part of same
    %mesh
    other_orig_elem = [];  other_orig_vert = [];
    for other_ind = other_inds
        other_orig_elem = [other_orig_elem; original_indices(other_ind).elem];
        other_orig_vert = [other_orig_vert; original_indices(other_ind).vert];
    end
    
    
    original_indices(ind).elem = original_indices(ind).elem + max(other_orig_elem);
    original_indices(ind).vert = original_indices(ind).vert + max(other_orig_vert);
    Mesh(ind).elems             = Mesh(ind).elems             + max(other_orig_vert);  %still using orig vert labels here, until renumber_Mesh below
end