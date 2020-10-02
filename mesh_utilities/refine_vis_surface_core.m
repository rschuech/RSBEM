function [Surface_vis] = refine_vis_surface_core(Mesh, n_refines, Surface_vis)

% for imesh = 1:length(Mesh)
%     Surface_vis_temp(imesh).nodes = Surface_vis(imesh).nodes;  %prevents parfor Coder complaint about variable classification
% end

Surface_vis_temp = Surface_vis;

parfor imesh = 1:length(Mesh)
    
    temp_surfun = [];
    
    if n_refines(imesh) == 0 %don't refine at all
        
        vis_elements = Mesh(imesh).elements(:,1:3); %remove midpoint nodes
        
        % have to do some work to get rid of now unused nodes but keep
        % numbering correct
        node_inds = unique(vis_elements(:));  %new node inds used in flat elements
        temp_elements = vis_elements;
        temp_nodes = NaN(length(node_inds),3);
        if isfield(Surface_vis_temp,'surfun')
        temp_surfun = NaN(length(node_inds),size(Surface_vis_temp(imesh).surfun,2));
        end
        
        c = 0;
        for vi = 1:length(node_inds)
            c = c + 1;
            temp_elements(temp_elements == node_inds(vi)) = c;
            temp_nodes(c,:) = Surface_vis_temp(imesh).nodes(node_inds(vi),:);
             if isfield(Surface_vis_temp,'surfun')
            temp_surfun(c,:) = Surface_vis_temp(imesh).surfun(node_inds(vi),:);
             end
        end
        
        
        Surface_vis(imesh).elements = temp_elements;
        Surface_vis(imesh).nodes = temp_nodes;
        Surface_vis(imesh).n_node = size(Surface_vis_temp(imesh).nodes,1); %recalc n_nodes
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