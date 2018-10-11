function [Mesh] = renumber_Mesh(Mesh)
%renumbers vert indices to start at 1 for each submesh, in terms of when
%they are referred to in each element
%leaves global indices and orig indices unchanged

%after running this, can refer to elements simply by their row index (trivial
%since nothing needs to really be renumbered), and verts can also be
%referred to by their row index (since below, we will renumber the element
%definitions to use this new vert numbering scheme)


for i = 1:length(Mesh)  %submeshes
    
    temp_elems = NaN(size(Mesh(i).elems)); %using a new copy of elems ensures we don't screw up revised vert inds we updated earlier in the loop
    for vert_i = 1:Mesh(i).n_vert %new vert indices always start at 1 and go to # verts in this submesh
        orig_ind = Mesh(i).indices.orig.vert(vert_i); %this is what is still referenced in the elements
        %the new vert index will simply be it's current row index, vert_i
        %do a replacement everywhere this vert appears in the elements
        temp_elems(Mesh(i).elems == orig_ind) = vert_i;
    end
    
    Mesh(i).elems = temp_elems;
    
    if any(isnan(Mesh(i).elems(:)))
        disp('Error:  renumbering of verts didn''t work')
        pause
    end
    
end