function [vert_range, elem_range] = bounds_from_global_inds(Mesh)
%outputs ranges of global indices for each submesh:
% each col of vert_range or elem_range is for a submesh,
%[starting index 1, starting index 2, ... ; ending index 1, ending index 2, ...]

vert_starts = NaN(1,length(Mesh));
vert_ends = vert_starts;  elem_starts = vert_starts;  elem_ends = vert_starts;

for i = 1:length(Mesh)  %submeshes
    %vert_starts(i) = Mesh(i).indices.global.vert(find(~isnan(Mesh(i).indices.global.vert),1,'first')); %could just use min() I guess but this should be faster?
    vert_starts(i) = Mesh(i).indices.glob.unq_bounds.vert(1);
    %global inds had better be monotonically increasing, but there may
    %conceivably be some NaNs at the beginning for skipped elements from
    %prior submeshes
    %vert_ends(i) = Mesh(i).indices.global.vert(find(~isnan(Mesh(i).indices.global.vert),1,'last')); %could just use max() I guess but this should be faster?
    vert_ends(i) = Mesh(i).indices.glob.unq_bounds.vert(2);
    %same as above, except now we want the last non-NaN entry since it must
    %be the largest
    %elem_starts(i) = Mesh(i).indices.global.elem(find(~isnan(Mesh(i).indices.global.elem),1,'first'));
    elem_starts(i) = Mesh(i).indices.glob.unq_bounds.elem(1);
   % elem_ends(i) = Mesh(i).indices.global.elem(find(~isnan(Mesh(i).indices.global.elem),1,'last'));
    elem_ends(i) = Mesh(i).indices.glob.unq_bounds.elem(2);
end


vert_range = [vert_starts; vert_ends];
elem_range = [elem_starts; elem_ends];
