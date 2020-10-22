y0 = [ Network0.l; reshape(Network0.nodes',Network0.n_nodes*3,1) ; Mesh0(2).refpoints; zeros(4,1)];

Mesh = Mesh0; Network = Network0; y = y0;

dt = 0.001;
t = 0:dt:1;

for c = 2:length(t)
    
    
    figure(125);
    plot3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),'ko','MarkerFaceColor','k');
    % text(Network.nodes(:,1)+0.2,Network.nodes(:,2),Network.nodes(:,3),cellfun(@num2str,num2cell(1:Network.n_nodes),'UniformOutput',false),'FontSize',12);
    hold on
    plot3([Network.nodes(Network.links(:,1),1) Network.nodes(Network.links(:,2),1)]',...
        [Network.nodes(Network.links(:,1),2) Network.nodes(Network.links(:,2),2)]',...
        [Network.nodes(Network.links(:,1),3) Network.nodes(Network.links(:,2),3)]',...
        'r--','LineWidth',2);
    
    % quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.u(:,1),Network.u(:,2),Network.u(:,3),'b-','LineWidth',2.5);
    % quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.g(:,1),Network.g(:,2),Network.g(:,3),'b--','LineWidth',2.5);
    
    [s,e] = plot_mesh(Mesh,1);  set(s,'facealpha',0.6);
    if c == 2
        
        xlabel('x'); ylabel('y'); zlabel('z');
    end
    
    hold off
    light
    
    axis equal
    % view(90,0);
    view(50, 0);
%      pause
drawnow
    pause
    
    
    
    
    
    [ A_temp, A_force_recycle, A_torque_recycle, RHS, A_motor_torque0] = matrix_assembly_mex_wrapper(Mesh,Network, matrix_props,index_mapping,mesh_node_parameters,assembly_input);
    
    [solution] = matrix_solve(assembly_input, matrix_props, A_temp, RHS);
    
    der_local = solution( matrix_props.n_collocation*3 + 1  : end);  % need to fix for free slip?
    
    
    [u_network_nodes] = field_velocity(Mesh,Network, Network.nodes, solution, matrix_props,index_mapping,mesh_node_parameters,assembly_input);
    
    
    
    dldt =  Network.E .* Network.l_0 ./ Network.eta .* ( r ./ Network.l - 1);
    
    
    dydt = [ dldt; reshape(u_network_nodes', Network.n_nodes*3,1); der_local ];
    
%     y(1:Network.n_nodes*3) = 
%     y = y + [dldt*dt; reshape(u_network_nodes', Network.n_nodes*3,1);
y(1:(Network.n_links + Network.n_nodes*3)) = y(1:(Network.n_links + Network.n_nodes*3)) + dydt(1:(Network.n_links + Network.n_nodes*3)) * dt;
y((Network.n_links + Network.n_nodes*3+1):end) = dydt((Network.n_links + Network.n_nodes*3+1):end) * dt;

%     y = y + dydt * dt;
% y(end-2:end-1) = 0;
    %% Update Mesh and Network based on y
%     [Mesh, Network, mesh_node_parameters, y_swimmer, r, A_1] = update_Mesh_Network(Mesh, Network, t, y, mesh_node_parameters,index_mapping, assembly_input);
    
der_local(1:3) = 0;

    Mesh2(1).nodes = Mesh(1).nodes + dt * (repmat(der_local(1:3)',Mesh(1).n_nodes,1) + cross( repmat(der_local(4:6)',Mesh(1).n_nodes,1) , Mesh(1).nodes - Mesh(2).refpoints' ) );
    Network.nodes = Network.nodes + dt * u_network_nodes;
    
    
%     Mesh3 = shiftMesh(Mesh, der_local(1:3)*dt);
    
    [Mesh3, mesh_node_parameters, A_1] = rotateMesh(Mesh,Mesh(2).refpoints(:,1),mesh_node_parameters,index_mapping, ...
    1, der_local([4 5 6])*dt);
    
    
end