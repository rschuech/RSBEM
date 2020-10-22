


switch input.Tail.motorBC
    case "torque"
        Mesh(1).u = zeros(Mesh(1).n_nodes,3);  Mesh(2).u = zeros(Mesh(2).n_nodes,3);
        assembly_input.Tail.motor_torque = input.Tail.motor_torque;
    case "freq"
        % need to set zero u on body, and relative u due to known freq on tail....
end

y0 = [ Network0.l; reshape(Network0.nodes',Network0.n_nodes*3,1) ; 0*Mesh0(1).refpoints; zeros(4,1)];

% Mesh0 = Mesh;  Network0 = Network;

% [dydt] = derivatives(0,y0, Mesh, Network, mesh_node_parameters, index_mapping, matrix_props, assembly_input);

mesh_node_parameters.normals_avg = NaN(size(mesh_node_parameters.normals_avg));
mesh_node_parameters.tangents_avg = NaN(size(mesh_node_parameters.tangents_avg));

fun = @(t,y) derivatives(t,y, Mesh0, Network0, mesh_node_parameters, index_mapping, matrix_props, assembly_input);


% tspan = [0 0.15];
tspan = [0 0.15];  % tail rotation is about 0.002 and 0.02 is about 10 tail rotations
% ops = odeset('OutputFcn',@odetpbar);
options = odeset('Refine',8,'RelTol',1E-3,'AbsTol',1E-3,'OutputFcn',@odetpbar);
[sol] = ode113(fun, tspan, y0,options);

% sol = odextend(sol, fun, 0.02);

tfinals = 0.2:0.05:3;
for tfinal = tfinals
    sol = odextend(sol, fun, tfinal);
end


%%
save_vid = true;
if save_vid
      try, close(vidh); end
    vidh = VideoWriter('C:\Hull\swimmer network videos\network swimmer test E 1 eta 10','MPEG-4');
    vidh.FrameRate = 30; %50;
    vidh.Quality = 95; %1-100
    resolution = '-r0';
    %     resolution = '-r180';  % doesn't seem to look much different, but
    %     doubles file size
  
    open(vidh);
end

tfinal = sol.x(end);
% tfinal = 0.003;

T = linspace(tspan(1),tfinal,2000); % 200

for t = T(1:end)
    
    y = deval(sol, t);
    
    [Mesh, Network, mesh_node_parameters, y_swimmer] = update_Mesh_Network(Mesh0, Network0,t, y, mesh_node_parameters,index_mapping, assembly_input);
    
    
    hfig = figure(133);
    n1 = 64;
    plot3(Network.nodes(1:n1,1),Network.nodes(1:n1,2),Network.nodes(1:n1,3),'ko','MarkerFaceColor','k','markersize',5);
    hold on
    plot3(Network.nodes(n1+1:end,1),Network.nodes(n1+1:end,2),Network.nodes(n1+1:end,3),'bo','MarkerFaceColor','b','markersize',3);
    % text(Network.nodes(:,1)+0.2,Network.nodes(:,2),Network.nodes(:,3),cellfun(@num2str,num2cell(1:Network.n_nodes),'UniformOutput',false),'FontSize',12);
    hold on
    plot3([Network.nodes(Network.links(:,1),1) Network.nodes(Network.links(:,2),1)]',...
        [Network.nodes(Network.links(:,1),2) Network.nodes(Network.links(:,2),2)]',...
        [Network.nodes(Network.links(:,1),3) Network.nodes(Network.links(:,2),3)]',...
        'r--','LineWidth',2);
    
    % quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.u(:,1),Network.u(:,2),Network.u(:,3),'b-','LineWidth',2.5);
    % quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.g(:,1),Network.g(:,2),Network.g(:,3),'b--','LineWidth',2.5);
    
    [s,e] = plot_mesh(Mesh,1);
    if t == tspan(1)
        
        xlabel('x'); ylabel('y'); zlabel('z');
    end
    
    hold off
    light
    
    axis equal
    view(0,90);
%     view(-55, 15);
view([45 15]);
    xlim([-2.5 7]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
    
    title("t = " + num2str(t))
    drawnow
%     pause
    
    if save_vid
        orig_mode = get(hfig, 'PaperPositionMode');
        
        set(hfig, 'PaperPositionMode', 'auto');
        cdata = print('-RGBImage','-opengl',resolution);
        % Restore figure to original state
        set(hfig, 'PaperPositionMode', orig_mode); % end
        
        %For the "OpenGL" renderer you can write a similar code. This technique will not work for the "painters" renderer.
        
        %Next, replace the use of GETFRAME from your code with IM2FRAME as follows:
        
        F = im2frame(cdata);
        
        writeVideo(vidh,F);
    end
    
    
end

if save_vid
    close(vidh);
end