function [Mesh, Network, mesh_node_parameters, y_swimmer, r] = update_Mesh_Network(Mesh, Network, t, y, mesh_node_parameters,index_mapping, assembly_input)

%% Mesh

y_swimmer = y( Network.n_links + Network.n_nodes*3 + 1  : end) ;

if length(y_swimmer) == 6 %constant rotation rate condition, add y(7) internally even though we're not solving for it
    y_swimmer =  [y_swimmer; t * assembly_input.Tail.motor_freq];  %assumes initial angle = 0.  This is a monotonically increasing phase angle, not yet wrapped to 0 - 2*pi
else
    %     y_swimmer = y_in;
    
end
y_swimmer(7) = mod(y_swimmer(7) + 2*pi, 2*pi);  % shift periodic tail phase angle back to 0 - 2*pi

% [Mesh, mesh_node_parameters] = rotateTail(Mesh, mesh_node_parameters,index_mapping,assembly_input.Tail.submesh_index, y_swimmer(7));
%
% shift tail refpoint to origin, rotate tail around tail axis, shift back to original tail refpoint
[Mesh, mesh_node_parameters] = rotateMesh(Mesh,Mesh(assembly_input.Tail.submesh_index).refpoints(:,1),mesh_node_parameters,index_mapping, ...
    assembly_input.Tail.submesh_index, y_swimmer(7),Mesh(assembly_input.Tail.submesh_index).orientation(:,1));

% shift body refpoint to origin, rotate entire mesh around x, y, z, shift back
[Mesh, mesh_node_parameters] = rotateMesh(Mesh,Mesh(2).refpoints(:,1),mesh_node_parameters,index_mapping, ...
    1:length(Mesh), y_swimmer(4:6));

% shift entire mesh in x, y, z.  assuming here that y_swimmer(1:3) represents current location of body refpoint
% Mesh = shiftMesh(  Mesh, ( y_swimmer(1:3) - Mesh(2).refpoints(:,1) )  );
% Mesh = shiftMesh(  Mesh,  y_swimmer(1:3) - Mesh(1).refpoints(:,1)    );
Mesh = shiftMesh(  Mesh,  y_swimmer(1:3)   );

% [Mesh,mesh_node_parameters] = move_Mesh(Mesh, mesh_node_parameters, index_mapping, y_swimmer(1:6)   );




%% Network
Network.l = y(1:Network.n_links);
Network.nodes = reshape( y(Network.n_links + 1 : Network.n_links + Network.n_nodes*3) , 3,[])';

d = Network.nodes(Network.links(:,2),:)  - Network.nodes(Network.links(:,1),:) ;
r = sqrt(sum( d.^2 , 2 ) );


% mindists = NaN(Network.n_links,1);
% tic

for n = 1:Network.n_links
    
    if Network.E(n) ~= 0
        
        
        edge = createEdge3d(Network.nodes( Network.links(n,1) , :), Network.nodes( Network.links(n,2) , :));
        
%         [dists] = distancePointEdge3d(Mesh(1).nodes, edge); % just body
        
          [dists] = distancePointEdge3d(vertcat(Mesh(assembly_input.link_breakage_submeshes).nodes), edge); % body and tail
        
        % mindists(n) = min(dist);
        
        if min(dists) < assembly_input.link_breakage_distance
            Network.E(n) = 0;
        end
        
    end
    
end
% toc


f_s = Network.l_0.^2 .* Network.E .* (r./Network.l - 1) .* d ./ r;

g = NaN(Network.n_nodes,3); % net Maxwell element viscoelastic force on each network node due to all attached links

parfor i = 1:Network.n_nodes
    g(i,:) = sum(   Network.link_members{i}(2,:)' .*  f_s( Network.link_members{i}(1,:) , :) , 1);
    
end

Network.g = g;  %clear g;


Network.PE = sum( Network.l_0.^3 .* Network.E .* (r./Network.l - 1).^2 );
