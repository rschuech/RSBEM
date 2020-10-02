function [Mesh] = shiftMesh(Mesh,shift)

%translate mesh in space by [shift(1) shift(2) shift(3)]

shift = shift(:)';

for i = 1:length(Mesh) %do for all submeshes in the input
    
    Mesh(i).nodes = Mesh(i).nodes + repmat(shift,Mesh(i).n_nodes,1);
    
    
    %     Mesh(i).centroid = Mesh(i).centroid + repmat(shift',1,Mesh(i).n_elem);
    Mesh(i).centroid = Mesh(i).centroid + repmat(shift',1,size(Mesh(i).centroid,2));
    if isfield(Mesh,'refpoints')
        Mesh(i).refpoints = Mesh(i).refpoints +  repmat(shift',1,size(Mesh(i).refpoints,2));
    end
  
    
    if isfield(Mesh,'Centroid')
        Mesh(i).Centroid = Mesh(i).Centroid + shift';
    end
    
end