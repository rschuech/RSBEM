function [Surface_vis] = refine_vis_surface_core(Mesh, n_refines, Surface_vis)



parfor imesh = 1:length(Mesh)
    
    
    
    if n_refines(imesh) == 0 %don't refine at all
        
        Surface_vis(imesh).elems = Mesh(imesh).elems(:,1:3); %remove midpoint nodes
        
        % have to do some work to get rid of now unused verts but keep
        % numbering correct
        vert_inds = unique(Surface_vis(imesh).elems(:));  %new vert inds used in flat elems
        temp_elems = Surface_vis(imesh).elems;
        temp_verts = NaN(length(vert_inds),3);
        
        c = 0;
        for vi = 1:length(vert_inds)
            c = c + 1;
            temp_elems(temp_elems == vert_inds(vi)) = c;
            temp_verts(c,:) = Surface_vis(imesh).verts(vert_inds(vi),:);
        end
        
        
        Surface_vis(imesh).elems = temp_elems;
        Surface_vis(imesh).verts = temp_verts;
        Surface_vis(imesh).n_vert = size(Surface_vis(imesh).verts,1); %recalc n_verts

    else % refining at least once
        temp = refine_wrap(Mesh(imesh), imesh,  n_refines);
        fields = fieldnames(temp);
        for f = 1:length(fields) %need to do this to avoid complaints about output of refine_wrap not having all the fields of input Mesh
            Surface_vis(imesh).(fields{f}) = temp.(fields{f});
        end
        
    end
    
end