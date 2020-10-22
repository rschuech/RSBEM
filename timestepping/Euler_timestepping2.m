

y0 = [ Network0.l; reshape(Network0.nodes',Network0.n_nodes*3,1) ; 0* Mesh0(2).refpoints; zeros(4,1)];

Mesh = Mesh0; Network = Network0; y = y0;

fun = @(t,y) derivatives(t,y, Mesh, Network, mesh_node_parameters, index_mapping, matrix_props, assembly_input);


dt = 0.00001;
t = 0:dt:0.1;

for c = 2:length(t)
    
    
    
figure(130);
plot3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),'ko','MarkerFaceColor','k');
% text(Network.nodes(:,1)+0.2,Network.nodes(:,2),Network.nodes(:,3),cellfun(@num2str,num2cell(1:Network.n_nodes),'UniformOutput',false),'FontSize',12);
hold on
plot3([Network.nodes(Network.links(:,1),1) Network.nodes(Network.links(:,2),1)]',...
    [Network.nodes(Network.links(:,1),2) Network.nodes(Network.links(:,2),2)]',...
    [Network.nodes(Network.links(:,1),3) Network.nodes(Network.links(:,2),3)]',...
    'r--','LineWidth',2);

% quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.u(:,1),Network.u(:,2),Network.u(:,3),'b-','LineWidth',2.5);
% quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.g(:,1),Network.g(:,2),Network.g(:,3),'b--','LineWidth',2.5);

[s,e] = plot_mesh(Mesh,1);
% if t == tspan(1)
    
    xlabel('x'); ylabel('y'); zlabel('z');
% end

hold off
light

axis equal
% view(90,0);
view(-19, 7.6);
 drawnow
 pause
 
 
 
 
 
 
    y(:,c) = y(:,c-1) + dt * fun(NaN, y(:,c-1));
% end


% for c = 1:length(t)

[Mesh, Network, mesh_node_parameters, y_swimmer] = update_Mesh_Network(Mesh0, Network0,NaN, y(:,c), mesh_node_parameters,index_mapping, assembly_input);


 
end

