function [Repulsion] = calc_repulsive_forces(Mesh, Network, index_mapping, repulsion)
% tic


repulsion.mindist2 = [0 0; 0 0];  % should cause dist2mesh to skip over everything fast

[distance2mesh, x, xi_eta, mesh_index, element_index, is_inside] = dist2mesh(Network.nodes, Mesh, index_mapping, repulsion.mindist2);



% toc
too_close = (distance2mesh < repulsion.d) | is_inside;
% Repulsion.distance2mesh(~too_close) = 0;  Repulsion.x(~too_close,:) = NaN;  Repulsion.mesh_index(~too_close) = 0;

% c = 0;
active = true(Network.n_nodes,1); % active if at least one working link, passive if no active links (all E == 0)
% active_close = true(sum(too_close),1);
for i = 1:Network.n_nodes
    active(i) = any( Network.E(  Network.link_members{i}(1,:) ) ~= 0 ); % node is active if it has at least one link with nonzero stiffness
    %     if too_close(i)
    %         c = c + 1;
    % %         active_close(c) = active(i); %surely a better way of doing this but whatever
    %         % length(active_close) == sum(too_close)
    %     end
end



% we will still keep track of both active and passive nodes here
Repulsion.network_node_indices = find(too_close);
Repulsion.mesh_index = mesh_index(too_close);
Repulsion.distance2mesh = distance2mesh(too_close);
Repulsion.x = x(too_close,:);


Repulsion.is_inside = is_inside(too_close);

% however, we will only impose repulsive forces on active nodes with links
% that are too close
too_close_active = too_close & active;

vec = Network.nodes(too_close_active,:) - x(too_close_active,:);
dir = vec ./ sqrt(  sum( ( vec ).^2 , 2 ) );
Repulsion.F = zeros(Network.n_nodes,3);
% repulsion forces on network nodes
Repulsion.F(too_close_active, :) = ...
    repulsion.g * ( exp(-distance2mesh(too_close_active)/repulsion.d) - exp(-1) ) ./ (1 - exp(-distance2mesh(too_close_active)/repulsion.d)) .* dir .* (1 - 2*is_inside(too_close_active));
% if is_inside = 0, then we multiply by 1; if is_inside = 1, then we multiply by -1, meaning we actually have a force pulling the offending network node toward the
% mesh, in the outward direction in hopes it will exit soon

