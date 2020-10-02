function [Mesh] = global_inds(Mesh)

%orig inds is what came directly from the .dat file.  however, orig inds may be repeated among submeshes, and this is a
%problem when assembling A matrix.  So make global index lists, which
%monotonically increase from submesh to submesh, and skip over elems or
%nodes that have already appeared in a previous submesh.  this should make
%assembling A easier....

concat.elem = [];  %concatenated elem inds
concat.node = [];

for i = 1:length(Mesh)  %submeshes
    
    names = {'elem','node'};
    for n = 1:length(names)
        name = names{n};
        %convention here is that elems, nodes are included in the first submesh in which
        %they occur (so if a node belongs to body and tail, it will only be
        %listed in body since it's the first submesh)
        switch name
            case 'elem'
                Mesh(i).indices.glob.elem = NaN(Mesh(i).n_elems,1);
            case 'node'
                Mesh(i).indices.glob.node = NaN(Mesh(i).n_nodes,1);
        end
        
        if i > 1 %at least 2nd submesh, concatenate prior submesh orig indices onto big aggregated list
            concat.(name) = [concat.(name); Mesh(i-1).indices.orig.(name)];
        end
        
        %unq_bounds contains the start and end global indices of elems or nodes
        %first found in this submesh, as opposed to previous submeshes (it does
        %not mean they are *only* found in the current submesh)
        
        %global indices outside the range of unq_bounds must be from
        %previous meshes in which the same indices were found before.  So
        %if when going through global indices, an index comes up that is
        %outside the range of unq_bounds, it should be skipped if we're
        %looping over collocation points, since it was already accounted
        %for
        
        if i == 1 %first submesh, start global inds at 1
            Mesh(i).indices.glob.unq_bounds.(name)(1) = 1;
        else
            Mesh(i).indices.glob.unq_bounds.(name)(1) = Mesh(i-1).indices.glob.unq_bounds.(name)(2) + 1;  %start at last final global ind + 1
        end
        
        [is_repeated, ~] = ismember(Mesh(i).indices.orig.(name), concat.(name));
        Mesh(i).indices.glob.unq_bounds.(name)(2) = (Mesh(i).indices.glob.unq_bounds.(name)(1) + sum( ~is_repeated) - 1);
        
        Mesh(i).indices.glob.(name)( ~is_repeated) = Mesh(i).indices.glob.unq_bounds.(name)(1) : Mesh(i).indices.glob.unq_bounds.(name)(2);
      
        repeated_inds = find(is_repeated);
        for ri = 1:length(repeated_inds)
            for mi = 1:(i-1) %start looking through previous submeshes for this orig index
                match = find(Mesh(mi).indices.orig.(name) == Mesh(i).indices.orig.(name)(repeated_inds(ri))); %local index in Mesh(mi) that matches this orig index
                if ~isempty(match)
                    Mesh(i).indices.glob.(name)(repeated_inds(ri)) = Mesh(mi).indices.glob.(name)(match);
                    break
                end
            end
        end
        
    end  %elem and node
    
end %submeshes

