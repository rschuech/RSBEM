function handles = label_elements(Mesh, inds)

%only labels elements going with indices (inds) for clarity

n = 0;

for j = inds  %element indices to label
    n = n+1;
    coords = Mesh.verts(Mesh.elems(j,1:3),:);
    center = mean(coords,1);
    labels(n) =  text(center(1),center(2),center(3),num2str(j),'fontsize',14);
    
end
