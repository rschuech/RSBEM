function [index_mapping] = create_global2local_index_mapping(Mesh, index_mapping)

%orig inds is what came directly from the .dat file.  however, orig inds may be repeated among submeshes, and this is a
%problem when assembling A matrix.  So make global index lists, which
%monotonically increase from submesh to submesh, and skip over elems or
%nodes that have already appeared in a previous submesh.  this should make
%assembling A easier....

n_global_nodes = max( vertcat( index_mapping.global2local.node ) );

% index_mapping.global2local.parent_submesh = cell( n_global_nodes , 1);
% index_mapping.global2local.node = cell( n_global_nodes , 1);


for submesh_ind = 1:length(index_mapping.local2global)  %submeshes
    
    for local_node_ind = 1:Mesh(submesh_ind).n_nodes
        global_ind = index_mapping.local2global(submesh_ind).node(local_node_ind);
       
%         index_mapping.global2local.parent_submesh{global_ind}(end+1) = submesh_ind;
%         index_mapping.global2local.node{global_ind}(end+1) = local_node_ind;
        [element_inds, position_inds] = find( Mesh(submesh_ind).elements == local_node_ind );
%         index_mapping.global2local.elements{global_ind}{end+1} = [element_inds  position_inds];
        
        % global2local is a n_global_nodes cell array, each cell is a matrix with rows = submesh, local node ind, element ind, position in element
        index_mapping.global2local.node{global_ind}(end+1:end+size(element_inds,1),:) = [repmat([submesh_ind, local_node_ind],size(element_inds,1),1), element_inds, position_inds];
    end
end
        index_mapping.global_node2local_node
        index_mapping.local_node2global_node
        index_mapping.global_node2global_unknown_u
        index_mapping.global_unknown_u2global_node
        
        
     
        


% concat.elem = [];  %concatenated elem inds
% concat.node = [];
% 
% for i = 1:length(Mesh)  %submeshes
%     
%     names = {'elem','node'};
%     for n = 1:length(names)
%         name = names{n};
%         %convention here is that elems, nodes are included in the first submesh in which
%         %they occur (so if a node belongs to body and tail, it will only be
%         %listed in body since it's the first submesh)
%         switch name
%             case 'elem'
%                 index_mapping.local2global(i).elem = NaN(Mesh(i).n_elems,1);
%             case 'node'
%                 index_mapping.local2global(i).node = NaN(Mesh(i).n_nodes,1);
%         end
%         
% %         if i > 1 %at least 2nd submesh, concatenate prior submesh orig indices onto big aggregated list
% %             concat.(name) = [concat.(name); Mesh(i-1).indices.orig.(name)];
% %         end
%         
%         %unq_bounds contains the start and end global indices of elems or nodes
%         %first found in this submesh, as opposed to previous submeshes (it does
%         %not mean they are *only* found in the current submesh)
%         
%         %global indices outside the range of unq_bounds must be from
%         %previous meshes in which the same indices were found before.  So
%         %if when going through global indices, an index comes up that is
%         %outside the range of unq_bounds, it should be skipped if we're
%         %looping over collocation points, since it was already accounted
%         %for
%         
% %         if i == 1 %first submesh, start global inds at 1
% %             Mesh(i).indices.glob.unq_bounds.(name)(1) = 1;
% %         else
% %             Mesh(i).indices.glob.unq_bounds.(name)(1) = Mesh(i-1).indices.glob.unq_bounds.(name)(2) + 1;  %start at last final global ind + 1
% %         end
%         
% %         [is_repeated, ~] = ismember(Mesh(i).indices.orig.(name), concat.(name));
% %         Mesh(i).indices.glob.unq_bounds.(name)(2) = (Mesh(i).indices.glob.unq_bounds.(name)(1) + sum( ~is_repeated) - 1);
%       
% temp = vertcat(Mesh.elems); 
% if length(unique(temp)) ~= max(temp(:))
%     error('Original global mesh indices are non-consecutive.  Must remove gaps and renumber.');
%     % if there are some missing global inds, the total number of nodes is not simply the largest global index. there would probably be other problems too.
%     % Most logical to simply make sure global indices are consecutive.
% end
% 
%  index_mapping.local2global(i).node = original_indices(i).node;
%  index_mapping.local2global(i).elem = original_indices(i).elem;
% 
%         Mesh(i).indices.glob.(name)( ~is_repeated) = Mesh(i).indices.glob.unq_bounds.(name)(1) : Mesh(i).indices.glob.unq_bounds.(name)(2);
%       
%         repeated_inds = find(is_repeated);
%         for ri = 1:length(repeated_inds)
%             for mi = 1:(i-1) %start looking through previous submeshes for this orig index
%                 match = find(Mesh(mi).indices.orig.(name) == Mesh(i).indices.orig.(name)(repeated_inds(ri))); %local index in Mesh(mi) that matches this orig index
%                 if ~isempty(match)
%                     Mesh(i).indices.glob.(name)(repeated_inds(ri)) = Mesh(mi).indices.glob.(name)(match);
%                     break
%                 end
%             end
%         end
%         
%     end  %elem and node
%     
% end %submeshes
% 
