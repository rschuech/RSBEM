function [Surface_vis] = refine_vis_surface_core(Mesh, n_refines, Surface_vis)

% for imesh = 1:length(Mesh)
%     Surface_vis_temp(imesh).verts = Surface_vis(imesh).verts;  %prevents parfor Coder complaint about variable classification
% end

Surface_vis_temp = Surface_vis;

parfor imesh = 1:length(Mesh)
    
    temp_surfun = [];
    
    if n_refines(imesh) == 0 %don't refine at all
        
        vis_elems = Mesh(imesh).elems(:,1:3); %remove midpoint nodes
        
        % have to do some work to get rid of now unused verts but keep
        % numbering correct
        vert_inds = unique(vis_elems(:));  %new vert inds used in flat elems
        temp_elems = vis_elems;
        temp_verts = NaN(length(vert_inds),3);
        if isfield(Surface_vis_temp,'surfun')
        temp_surfun = NaN(length(vert_inds),size(Surface_vis_temp(imesh).surfun,2));
        end
        
        c = 0;
        for vi = 1:length(vert_inds)
            c = c + 1;
            temp_elems(temp_elems == vert_inds(vi)) = c;
            temp_verts(c,:) = Surface_vis_temp(imesh).verts(vert_inds(vi),:);
             if isfield(Surface_vis_temp,'surfun')
            temp_surfun(c,:) = Surface_vis_temp(imesh).surfun(vert_inds(vi),:);
             end
        end
        
        
        Surface_vis(imesh).elems = temp_elems;
        Surface_vis(imesh).verts = temp_verts;
        Surface_vis(imesh).n_vert = size(Surface_vis_temp(imesh).verts,1); %recalc n_verts
         if isfield(Surface_vis_temp,'surfun')
        Surface_vis(imesh).surfun = temp_surfun;
         end

    else % refining at least once
        temp = refine_wrap(Mesh(imesh), imesh,  n_refines, Surface_vis_temp(imesh));
        fields = fieldnames(temp);
        for f = 1:length(fields) %need to do this to avoid complaints about output of refine_wrap not having all the fields of input Mesh
            Surface_vis(imesh).(fields{f}) = temp.(fields{f});
        end
        
    end
    
end