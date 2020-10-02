function edges = get_edges(Mesh)

for m = 1:length(Mesh)
    
    edges(m).verts = zeros(0,2);  %list of unique edges defined by two vertices always given from small to large vert ind
    edges(m).elems = {};  %list of all the elements that share each edge (normally 2 triangles per edge unless something is majorly wrong!)
    
    for i = 1:Mesh(m).n_elem
        
        elem = sort( Mesh(m).elems(i,1:3) );  %corner verts only
        
        edges_temp = [elem(1) elem(2); elem(1) elem(3); elem(2) elem(3)]; %always small then big index
        
        [ispresent, inds] = ismember(edges_temp, edges(m).verts,'rows');
        
        edges(m).verts( end+1:end+sum(~ispresent) ,:) = edges_temp(~ispresent,:);
        
        if any(ispresent)
            for j = inds(ispresent)'
                edges(m).elems{j,1}(end+1) = i;
            end
        end
        
        for j = 1:sum(~ispresent)
            edges(m).elems{end+1,1} = i;
        end
        
    end
    
end


%%
% figure
%
% for i = 1:size(edges.verts,1)
%     vert1 = Mesh.verts(edges.verts(i,1),:);
%     vert2 = Mesh.verts(edges.verts(i,2),:);
%     plot3( [vert1(1) vert2(1)],[vert1(2) vert2(2)],[vert1(3) vert2(3)],'-r');
%     hold on
% end
% axis equal
% grid
