function [combined] = combine_Meshes(Mesh,Inds)
% combines submeshes, getting rid of duplicated verts and elems and just
% outputting a list of global verts (always the same) and a list of
% elements for only Inds submeshes, so operations can be performed on this
% combined smaller mesh

global_verts = [];
for i = 1:length(Mesh)
    inds = find( Mesh(i).indices.glob.vert >= Mesh(i).indices.glob.unq_bounds.vert(1) &   Mesh(i).indices.glob.vert <= Mesh(i).indices.glob.unq_bounds.vert(2)); %local inds of "global" verts belonging to this mesh
    global_verts = [global_verts; Mesh(i).verts(inds,:)];
end



elemlist = [];
for i = Inds
    temp = Mesh(i).elems;
    for j = 1:Mesh(i).n_vert
        temp(Mesh(i).elems == j) = Mesh(i).indices.glob.vert(j);
    end
    elemlist = [elemlist; temp];
end


combined.verts = global_verts;
combined.elems = elemlist;