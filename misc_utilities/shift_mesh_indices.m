function [Mesh] = shift_mesh_indices(names,Mesh)

ind = find(strcmp('Tail',names));  %Mesh ind going with tail, if we loaded tail
if ~isempty(ind)
    other_inds = setdiff(1:length(names), ind);  % inds for all other submeshes
    %tail indices might overlap body/transverse/wingtip since they're not part of same
    %mesh
    other_orig_elem = [];  other_orig_vert = [];
    for other_ind = other_inds
        other_orig_elem = [other_orig_elem; Mesh(other_ind).indices.orig.elem];
        other_orig_vert = [other_orig_vert; Mesh(other_ind).indices.orig.vert];
    end
    
    
    Mesh(ind).indices.orig.elem = Mesh(ind).indices.orig.elem + max(other_orig_elem);
    Mesh(ind).indices.orig.vert = Mesh(ind).indices.orig.vert + max(other_orig_vert);
    Mesh(ind).elems             = Mesh(ind).elems             + max(other_orig_vert);  %still using orig vert labels here, until renumber_Mesh below
end