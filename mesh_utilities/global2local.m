function [i_mesh, local_index] = global2local(global_index, Mesh, index_range, type)

%convert from global to local inds for verts or elems
i_mesh = NaN(size(global_index));  local_index = i_mesh;

for i = 1:numel(global_index)
    
    temp =  find( global_index(i) >= index_range(1,:) & global_index(i) <= index_range(2,:) );
    if length(temp) > 1
        disp('Error:  "unique" global index found in more than one submesh')
        temp = NaN;
        pause
    end
    
    i_mesh(i) = temp(1); %appease Coder
        
    local_index(i) = find(Mesh(i_mesh(i)).indices.glob.(type) == global_index(i));
    
    
end
