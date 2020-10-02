function [Edge_vis] = refine_vis_edges(Mesh, sidepts)

%creates edge data for plotting curved triangles by creating many-sided polygonal non-planar elements that closely approximate curved triangles
% sidepts = # evenly sampled points along each side of each T6 triangle,
% either scalar that applies to each submesh or vector of values
% correponding to each submesh

nplaces = 9;  %for rounding away roundoff errors

if isscalar(sidepts)
    sidepts = repmat(sidepts,length(Mesh),1);
end


Edge_vis = struct('nodes',repmat({[]},length(Mesh),1), 'elements',repmat({[]},length(Mesh),1));
%     'n_vert',repmat({NaN},length(Mesh),1),'n_elements',repmat({NaN},length(Mesh),1) );  %initialize as input, all fields will be copied


parfor i_mesh = 1:length(Mesh)
    
    
    ii = 0;
    for elem_i = 1:Mesh(i_mesh).n_elements
        start = ii + 1;
        
        vertinds = Mesh(i_mesh).elements(elem_i,:);
        subnodes = Mesh(i_mesh).nodes(vertinds,:);
        shape_parameters = Mesh(i_mesh).shape_parameters(elem_i,:);  %[alpha beta gamma]
        
        
        %bottom side
        eta = 0;
        for Xi = linspace(0,1,sidepts(i_mesh))
            ii = ii+1;
            [Edge_vis(i_mesh).nodes(:,ii)] = T6interp(subnodes,Xi,eta,shape_parameters);
        end
        
        for Xi = linspace(1,0,sidepts(i_mesh))
            ii = ii+1;
            eta = -Xi + 1; %diagonal side
            [Edge_vis(i_mesh).nodes(:,ii)] = T6interp(subnodes,Xi,eta,shape_parameters);
        end
        
        xi = 0;
        for Eta = linspace(1,0,sidepts(i_mesh))
            ii = ii+1;
            [Edge_vis(i_mesh).nodes(:,ii)] = T6interp(subnodes,xi,Eta,shape_parameters);
        end
        
        last = ii;
        
        Edge_vis(i_mesh).elements(elem_i,:) = start:last;
        
    end
    
    Edge_vis(i_mesh).nodes = Edge_vis(i_mesh).nodes';
    
    % remove duplicated and shared nodes
    [un,~,ic] = unique(roundn(Edge_vis(i_mesh).nodes,-nplaces),'rows');
    %un = A(ia,:) and A = un(ic,:)
    
    temp = Edge_vis(i_mesh).elements;
    for i = 1:size(Edge_vis(i_mesh).elements,1)
        temp(i,:) = ic(Edge_vis(i_mesh).elements(i,:));
    end
    
    Edge_vis(i_mesh).nodes = un;
    Edge_vis(i_mesh).elements = temp;
    
end



for i_mesh = 1:length(Mesh)
    
    Edge_vis(i_mesh).n_vert = size(Edge_vis(i_mesh).nodes,1);
    Edge_vis(i_mesh).n_elements = size(Edge_vis(i_mesh).elements,1);
    
    fields = {'refpoints','orientation','name'};
    for f = 1:length(fields)
        if isfield(Mesh(i_mesh),fields{f})
            Edge_vis(i_mesh).(fields{f}) = Mesh(i_mesh).(fields{f});
        end
    end
    
end