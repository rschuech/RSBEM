function [corner_elems, flat_verts] = flatten_mesh(Mesh)

 % make new Faces and Vertices without midpoint nodes since these verts make
% point2trimesh angry
corner_elems = Mesh.elems(:,1:3);  % excludes midpoint nodes on curved triangles
flat_verts = [];
for ci = 1:numel(corner_elems)
    elem = corner_elems(ci);
    vert = Mesh.verts(elem,:);
    if ~isempty(flat_verts)
        [~, ind] = ismember(vert,flat_verts,'rows');
    else
        ind = 0;
    end

    if ~ind
        flat_verts(end+1,:) = vert;
        corner_elems(ci) = size(flat_verts,1);
    else
        corner_elems(ci) = ind;
    end
end