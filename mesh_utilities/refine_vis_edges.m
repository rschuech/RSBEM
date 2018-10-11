function [Edge_vis] = refine_vis_edges(Mesh, sidepts)

%creates edge data for plotting curved triangles by creating many-sided polygonal non-planar elements that closely approximate curved triangles
% sidepts = # evenly sampled points along each side of each T6 triangle,
% either scalar that applies to each submesh or vector of values
% correponding to each submesh

nplaces = 9;  %for rounding away roundoff errors

if isscalar(sidepts)
    sidepts = repmat(sidepts,length(Mesh),1);
end


Edge_vis = struct('verts',repmat({[]},length(Mesh),1), 'elems',repmat({[]},length(Mesh),1));
%     'n_vert',repmat({NaN},length(Mesh),1),'n_elem',repmat({NaN},length(Mesh),1) );  %initialize as input, all fields will be copied


parfor i_mesh = 1:length(Mesh)
    
    
    i = 0;
    for elem_i = 1:Mesh(i_mesh).n_elem
        start = i + 1;
        
        vertinds = Mesh(i_mesh).elems(elem_i,:);
        subverts = Mesh(i_mesh).verts(vertinds,:);
        shape_parameters = Mesh(i_mesh).elem_params(:,elem_i);  %[alpha beta gamma]
        
        
        %bottom side
        eta = 0;
        for xi = linspace(0,1,sidepts(i_mesh))
            i = i+1;
            [Edge_vis(i_mesh).verts(:,i)] = T6interp(subverts,xi,eta,shape_parameters);
        end
        
        for xi = linspace(1,0,sidepts(i_mesh))
            i = i+1;
            eta = -xi + 1; %diagonal side
            [Edge_vis(i_mesh).verts(:,i)] = T6interp(subverts,xi,eta,shape_parameters);
        end
        
        xi = 0;
        for eta = linspace(1,0,sidepts(i_mesh))
            i = i+1;
            [Edge_vis(i_mesh).verts(:,i)] = T6interp(subverts,xi,eta,shape_parameters);
        end
        
        last = i;
        
        Edge_vis(i_mesh).elems(elem_i,:) = start:last;
        
    end
    
    Edge_vis(i_mesh).verts = Edge_vis(i_mesh).verts';
    
    % remove duplicated and shared verts
    [un,~,ic] = unique(roundn(Edge_vis(i_mesh).verts,-nplaces),'rows');
    %un = A(ia,:) and A = un(ic,:)
    
    temp = Edge_vis(i_mesh).elems;
    for i = 1:size(Edge_vis(i_mesh).elems,1)
        temp(i,:) = ic(Edge_vis(i_mesh).elems(i,:));
    end
    
    Edge_vis(i_mesh).verts = un;
    Edge_vis(i_mesh).elems = temp;
    
end



for i_mesh = 1:length(Mesh)
    
    Edge_vis(i_mesh).n_vert = size(Edge_vis(i_mesh).verts,1);
    Edge_vis(i_mesh).n_elem = size(Edge_vis(i_mesh).elems,1);
    
    fields = {'refpoints','orientation','name'};
    for f = 1:length(fields)
        if isfield(Mesh(i_mesh),fields{f})
            Edge_vis(i_mesh).(fields{f}) = Mesh(i_mesh).(fields{f});
        end
    end
    
end