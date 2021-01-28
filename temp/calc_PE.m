PE = NaN(length(stored_output.time),1);
% is_inside = false(length(stored_output.time),Network0.n_nodes);
% distance2mesh = sparse(length(stored_output.time),Network0.n_nodes);
% submesh_index = sparse(length(stored_output.time),Network0.n_nodes);

t_repulsion = [];
is_inside = {};
distance2mesh = {};
submesh_index = {};
network_indices = {};
stored_output.repulsion.total_force_mag = zeros(length(stored_output.time),1);
stored_output.repulsion.total_torque_mag = zeros(length(stored_output.time),1);
stored_output.repulsion.any_inside = false(length(stored_output.time),1);
any_flagged = false(length(stored_output.time),1);

mesh_node_parameters0 = mesh_node_parameters;

outer_timer = tic;
for c = 1:length(stored_output.time)
    if rem(c,50) == 0
        c/length(stored_output.time)
    end
    t = stored_output.time(c);
    
    y = deval(sol, t);
    
    Network = Network0;
    
    Network.E(stored_output.link_breakage_time <= t) = 0;
    
    Network.l = y(1:Network0.n_links);
    Network.nodes = reshape( y(Network0.n_links + 1 : Network0.n_links + Network0.n_nodes*3) , 3,[])';
    
    
    %     [Mesh, Network, ~, y_swimmer] = update_Mesh_Network(Mesh0, Network,t, y, mesh_node_parameters,index_mapping, assembly_input);
    
    d = Network.nodes(Network.links(:,2),:)  - Network.nodes(Network.links(:,1),:) ;
    r = sqrt(sum( d.^2 , 2 ) );
    %     f_s = Network.l_0.^2 .* Network.E .* (r./Network.l - 1) .* d ./ r;
    %     f = sqrt(sum(f_s.^2,2))  .*  sign( r./Network.l - 1);
    PE(c) = sum( Network.l_0.^3 .* Network.E .* (r./Network.l - 1).^2 );
    
    
    
    
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
    [Mesh, mesh_node_parameters] = rotateMesh(Mesh0,Mesh0(assembly_input.Tail.submesh_index).refpoints(:,1),mesh_node_parameters0,index_mapping, ...
        assembly_input.Tail.submesh_index, y_swimmer(7),Mesh0(assembly_input.Tail.submesh_index).orientation(:,1));
    
    % shift body refpoint to origin, rotate entire mesh around x, y, z, shift back
    [Mesh, mesh_node_parameters] = rotateMesh(Mesh,Mesh(2).refpoints(:,1),mesh_node_parameters0,index_mapping, ...
        1:length(Mesh), y_swimmer(4:6));
    
    % shift entire mesh in x, y, z.  assuming here that y_swimmer(1:3) represents current location of body refpoint
    % Mesh = shiftMesh(  Mesh, ( y_swimmer(1:3) - Mesh(2).refpoints(:,1) )  );
    % Mesh = shiftMesh(  Mesh,  y_swimmer(1:3) - Mesh(1).refpoints(:,1)    );
    Mesh = shiftMesh(  Mesh,  y_swimmer(1:3)   );
    
    
    inner_timer  = tic;
    [Repulsion] = calc_repulsive_forces(Mesh, Network, index_mapping, repulsion);
    if rem(c,10) == 0
        toc(inner_timer);
    end
    
    
    
    if ~isempty(Repulsion.distance2mesh)
        t_repulsion(end+1,1) = t;
        any_flagged(c) = true;
        is_inside{end+1,1} = Repulsion.is_inside;
        distance2mesh{end+1,1} = Repulsion.distance2mesh;
        submesh_index{end+1,1} = Repulsion.mesh_index;
        network_indices{end+1,1} = Repulsion.network_node_indices;
        
        stored_output.repulsion.total_force_mag(c,1) = sqrt(sum( (sum(Repulsion.F(Repulsion.network_node_indices,:), 1)).^2 ))  ;
        [refpoint] = get_rotational_references(Mesh, assembly_input);
        
        stored_output.repulsion.total_torque_mag(c,1) = sqrt(sum(  sum( cross( (Repulsion.x  - refpoint') , Repulsion.F(Repulsion.network_node_indices,:) ), 1).^2  ));
        
        
        
    end
    
end

toc(outer_timer)

stored_output.PE = PE;
% stored_output.repulsion.t = t_repulsion;
stored_output.repulsion.any_flagged = any_flagged;
stored_output.repulsion.network_node_indices = network_indices;
stored_output.repulsion.is_inside = is_inside;
stored_output.repulsion.any_inside(stored_output.repulsion.any_flagged) = cellfun(@any, stored_output.repulsion.is_inside);
stored_output.repulsion.distance2mesh = distance2mesh;
stored_output.repulsion.submesh_index = submesh_index;