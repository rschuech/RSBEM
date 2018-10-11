function [Mesh] = shiftMesh(Mesh,shift)

%translate mesh in space by [shift(1) shift(2) shift(3)]

shift = shift(:)';

for i = 1:length(Mesh) %do for all submeshes in the input
    
    Mesh(i).verts = Mesh(i).verts + repmat(shift,Mesh(i).n_vert,1);
    
    
    Mesh(i).centroid = Mesh(i).centroid + repmat(shift',1,Mesh(i).n_elem);
    if isfield(Mesh,'refpoints')
        Mesh(i).refpoints = Mesh(i).refpoints +  repmat(shift',1,size(Mesh(i).refpoints,2));
    end
    
    if isfield(Mesh,'Centroid')
        Mesh(i).Centroid = Mesh(i).Centroid + shift';
    end
    
end