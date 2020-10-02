function [Mesh, local_node2global_node] = shift_global_indices(Mesh, local_node2global_node, coincident_submeshes)

names = [Mesh.name];

% ind = find(strcmp('Tail',names));  %Mesh ind going with tail, if we loaded tail
% if ~isempty(ind)
%     other_inds = setdiff(1:length(names), ind);  % inds for all other submeshes
%     %tail indices might overlap body/transverse/wingtip since they're not part of same
%     %mesh
%     other_orig_elem = [];  other_orig_vert = [];
%     for other_ind = other_inds
%         other_orig_elem = [other_orig_elem; original_indices(other_ind).elem];
%         other_orig_vert = [other_orig_vert; original_indices(other_ind).vert];
%     end
%
%
%     original_indices(ind).elem = original_indices(ind).elem + max(other_orig_elem);
%     original_indices(ind).vert = original_indices(ind).vert + max(other_orig_vert);
%     Mesh(ind).elems             = Mesh(ind).elems             + max(other_orig_vert);  %still using orig vert labels here, until renumber_Mesh below
% end


for c = 1:length(coincident_submeshes)
  %  if iscell( coincident_submeshes{c} ) % aggregate all submeshes in this coincident group, regardless of geometric entity organization
        [~,coincident_inds{c}] = ismember([coincident_submeshes{c}{:}], names);  % submesh indices going with each group of coincident submeshes
        % can't properly test this until I generate some meshes like this!
%         disp('check results, this is untested!');
%         pause
%     else % no cells means that each submesh is a different geometric entity
%         [~,coincident_inds{c}] = ismember(coincident_submeshes{c}, names);  % submesh indices going with each group of coincident submeshes
%     end
    temp =   vertcat( local_node2global_node{coincident_inds{c}} );
    max_ind(c) = max(temp);  % largest original index across all coincident submeshes in this group
    min_ind(c) = min(temp);
end

[~,sorted_coincident_groups] = sort(max_ind); % indices for coincident mesh groups sorted by increasing max node index
% max_ind = max_ind(sorted_coincident_groups);
% min_ind = min_ind(sorted_coincident_groups);


for c = 1:length(coincident_submeshes) - 1 % each group of coincident submeshes; the first group is not altered (note c+1 indexing in loop)
    cgi1 = sorted_coincident_groups(c); % lower coincident group index
    cgi2 = sorted_coincident_groups(c+1); % upper coincident group index
    increment(c+1) =  max_ind(cgi1) - min_ind(cgi2) + 1;  % what to increment all submeshes in cgi2 by to eliminate overlap with submeshes in cgi1
    for i = coincident_inds{cgi2} % all submeshes in cgi2
        local_node2global_node{i} = local_node2global_node{i} + sum(increment(1:c+1)); % must increment all submeshes in cgi2 by sum of
        % prior increments to take care of overlap with all prior submeshes
%         Mesh(i).elements = Mesh(i).elements + sum(increment(1:c+1)); % must do same thing in list of elements where global node indices are referenced
% no, don't want to change local element or node indices, right?!
    end
end

% at this point, any possible overlap in global indices is eliminated.  Theoretically there may be missing global indices (even though there shouldn't be),
% so let's renumber them to make sure that global indices is a list of consecutive integers with no gaps

unique_global_inds = [0; unique( vertcat( local_node2global_node{:} ) )]; % unique global indices, possibly with gaps  insert 0 as fake first index to account for 
% (unlikely) possiblity of global indices not starting at 1

diffs = diff(unique_global_inds);
gap_inds = find(diffs > 1);
local_node2global_node0 = local_node2global_node;  % store copy of original
Mesh0 = struct('elements',{Mesh.elements}); % store copy of original Mesh elements

for i = 1:length(gap_inds)
    for j = 1:length(Mesh)
        shift_down = local_node2global_node0{j} >= unique_global_inds(gap_inds(i)); % global inds to shift downward to eliminate this gap
        local_node2global_node{j}(shift_down) = local_node2global_node{j}(shift_down) - (diffs(gap_inds(i)) - 1); % shift any indices after the gap down, eliminating this gap
        % need to do same thing to element lists
        % no we don't, we want to keep element lists using the original local node indices, right?
%         shift_down = Mesh0(j).elements >= unique_global_inds(gap_inds(i)); % linear indices since elems is an N x 6 array
%         Mesh(j).elements(shift_down) = Mesh(j).elements(shift_down) - (diffs(gap_inds(i)) - 1);
    end
end

% now we can be sure that there is no overlap between global indices of different groups of coincident meshes, and also that the list of global indices is
% consecutive, starting at 1