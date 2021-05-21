


% switch input.Tail.motorBC
%     case "torque"
%         Mesh(1).u = zeros(Mesh(1).n_nodes,3);  Mesh(2).u = zeros(Mesh(2).n_nodes,3);
%         assembly_input.Tail.motor_torque = input.Tail.motor_torque;
%     case "freq"
%         % need to set zero u on body, and relative u due to known freq on tail....
% end


Mesh0 = Mesh;  Network0 = Network;
y0 = [ Network0.l; reshape(Network0.nodes',Network0.n_nodes*3,1) ; 0*Mesh0(1).refpoints(:,1); zeros(4,1)];

% Mesh0 = Mesh;  Network0 = Network;

% [dydt] = derivatives(0,y0, Mesh, Network, mesh_node_parameters, index_mapping, matrix_props, assembly_input);

mesh_node_parameters.normals_avg = NaN(size(mesh_node_parameters.normals_avg));
mesh_node_parameters.tangents_avg = NaN(size(mesh_node_parameters.tangents_avg));

assembly_input.performance.verbose = false;
assembly_input.Tail.motor_torque = @(t)motor_torque(t,torque0,1000);

fun = @(t,y) derivatives(t,y, Mesh0, Network0, mesh_node_parameters, index_mapping, matrix_props, assembly_input);


% tspan = [0 0.15];
tspan = [0 0.01];  % tail rotation is about 0.002 and 0.02 is about 10 tail rotations
% ops = odeset('OutputFcn',@odetpbar);
options = odeset('Refine',1,'RelTol',5E-4,'AbsTol',5E-4,'OutputFcn',@odetpbar);
clear textprogressbar
clear matrix_assembly_mex
%%
clear sol
[sol(1)] = ode45(fun, tspan, y0,options);

[sol(2)] = ode113(fun, [0.15 0.25], sol.y(:,end),options);

% sol = odextend(sol, fun, 0.02);

clear sol
clear textprogressbar
[sol] = ode45(fun, tspan, y0,options);

clear textprogressbar
tfinals = 0.025:0.0025:3;
for tfinal = tfinals
    sol = odextend(sol, fun, tfinal);
    save('c:\Users\rudi\Desktop\RD\temp\shat.mat','sol');
end


%%
Mesh0 = Mesh;  Network0 = Network;
y0 = [ Network0.l; reshape(Network0.nodes',Network0.n_nodes*3,1) ; 0*Mesh0(1).refpoints(:,1); zeros(4,1)];
assembly_input.performance.verbose = false;
assembly_input.Tail.motor_torque = @(t)motor_torque(t,torque0,1000);
clear textprogressbar
clear matrix_assembly_mex
solver = @ode45;

% 'RelTol',5E-4,'AbsTol',5E-4,  originally
% defaults are 1E-3, 1E-6
% timestepping_tol.reltol = 1E-6;
% timestepping_tol.abstol = 1E-10;
% timestepping_tol.reltol = 5E-4;
% timestepping_tol.abstol = 5E-4;
timestepping_tol.reltol = 1E-3;
timestepping_tol.abstol = 1E-8;
% timestepping_tol.reltol = 5E-4;
% timestepping_tol.abstol = 5E-4;

[sol, stored_output, trespassers ] = timestepping_wrapper([], [],y0, 0.01, Mesh0, Network0, mesh_node_parameters, index_mapping, matrix_props, assembly_input,timestepping_tol, solver);

clear textprogressbar
% tfinals = 0.005:0.005:2;
tfinals = 0.02:0.01:2;
for tfinal = tfinals
    if sol.x(end) >= tfinal
        continue
    end
   [sol, stored_output , trespassers] = timestepping_wrapper([sol], [stored_output],y0, tfinal, Mesh0, Network0, mesh_node_parameters, index_mapping, matrix_props, assembly_input,timestepping_tol, solver);
% if any(trespassers.is_inside)
  if any(trespassers.distance2mesh(trespassers.is_inside) >= 0.006 )
      pause
      
%       Network_temp = Network0;
%        Network_temp.E( ~isnan(stored_output.link_breakage_time) ) = 0;  % permanently set E = 0 each time a link breaks
% the hell with it, the links don't matter for just figuring out where the
% network nodes are relative to the swimmer when the event occurs
      [Mesh, Network, mesh_node_parameters, y_swimmer, r] = update_Mesh_Network(Mesh0, Network0, sol.xe(end), sol.ye(:,end), mesh_node_parameters,index_mapping, assembly_input);
      
        [distance2mesh, x, xi_eta, mesh_index, element_index, is_inside] = dist2mesh(Network.nodes, Mesh, index_mapping, assembly_input.repulsion.mindist2);
        
%         Network.l = y(1:Network.n_links);
network_nodes = reshape( sol.ye(Network.n_links + 1 : Network.n_links + Network.n_nodes*3, end) , 3,[])';
        network_nodes( is_inside,: ) = x(is_inside, :);

       network_nodes = [ reshape(network_nodes', Network.n_nodes*3,1) ];
       
 y0_new = sol.ye(:,end);
 
       y0_new( Network.n_links + 1 : Network.n_links + Network.n_nodes*3  ) = network_nodes;
        
       
       
       
       % last value of sol should be the same as the event, exactly at the
       % threshold, so ideally want to keep this, but then sol definitely
       % won't line up with stored_output (which won't have output at the
       % exact threshold since the solver doesn't actually take a step at
       % the event) so perhaps for now just overwrite the event moment in
       % sol, so we generally have the solution at the step before the
       % event, and then the right side of the discontinuity when the nodes
       % are moved onto the mesh
       % last value of stored_output should be the overshoot of the
       % threshold but it should get overwritten by default when we do
       % odextend
      
  end
% end
   %     save('c:\Users\rudi\Desktop\RD\temp\broken.mat','sol');
end

%%
 is_inside = vertcat(stored_output.repulsion.is_inside{:});
 dists = vertcat(stored_output.repulsion.distance2mesh{:});
 dists_inside = dists(is_inside);
 figure;  histogram(dists_inside)
 title('network nodes inside swimmer')
 xlabel('distance to mesh')
 
%%
tic  
save_vid = false;
if save_vid
      try, close(vidh); end
%     vidh = VideoWriter('C:\Hull\swimmer network videos\longer tail','MPEG-4');
     vidh = VideoWriter('C:\Users\rudi\Desktop\RD\Tulane vids\longer tail floppier network3','Motion JPEG AVI'); 
    vidh.FrameRate = 50; %50; 30
    vidh.Quality = 95; %1-100
    resolution = '-r0';
    %     resolution = '-r180';  % doesn't seem to look much different, but
    %     doubles file size
  
    open(vidh);
end

tfinal = sol.x(end);
%  tfinal = 0.06;

% T = linspace(0,tfinal,5); % 600  
% T = sol(1).x( sol(1).x <= Inf);
T = 0:1E-3:tfinal;
T = 0:5E-2:tfinal;

for t = T(  T >= 0 & T <= Inf)  %  sol(end).x(end)
    
    if t < 0.15
        ind = 1;
    else
        ind = 1;
    end
    
    y = deval(sol(ind), t);
    
    [Mesh, Network, mesh_node_parameters, y_swimmer] = update_Mesh_Network(Mesh0, Network0,t, y, mesh_node_parameters,index_mapping, assembly_input);
    
    
%     dists = sqrt( sum((Network.nodes - Mesh(1).refpoints(:,1)').^2 , 2) ) - Metadata(1).geom.sphererad;
%     insides = find(dists < 0);
%     rel_dists = dists / Metadata(1).geom.sphererad;
    
    
    hfig = figure(142);
    n1 = size(Network.nodes,1);
    pl1 = plot3(Network.nodes(1:n1,1),Network.nodes(1:n1,2),Network.nodes(1:n1,3),'ko','MarkerFaceColor','k','markersize',5);
    hold on
    pl2 = plot3(Network.nodes(insides,1),Network.nodes(insides,2),Network.nodes(insides,3),'bo','MarkerFaceColor','b','markersize',10);
%     plot3(Network.nodes(n1+1:end,1),Network.nodes(n1+1:end,2),Network.nodes(n1+1:end,3),'bo','MarkerFaceColor','b','markersize',3);
    % text(Network.nodes(:,1)+0.2,Network.nodes(:,2),Network.nodes(:,3),cellfun(@num2str,num2cell(1:Network.n_nodes),'UniformOutput',false),'FontSize',12);
    hold on
    plot3([Network.nodes(Network.links(:,1),1) Network.nodes(Network.links(:,2),1)]',...
        [Network.nodes(Network.links(:,1),2) Network.nodes(Network.links(:,2),2)]',...
        [Network.nodes(Network.links(:,1),3) Network.nodes(Network.links(:,2),3)]',...
        '--','LineWidth',0.01,'Color',repmat(0.75,1,3));
    
    % quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.u(:,1),Network.u(:,2),Network.u(:,3),'b-','LineWidth',2.5);
    % quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.g(:,1),Network.g(:,2),Network.g(:,3),'b--','LineWidth',2.5);
    
    [s,e] = plot_mesh(Mesh,2); s(1).FaceAlpha = 1;
    e(1).EdgeAlpha = 0.2;
    
    if t == 0
        
        xlabel('x'); ylabel('y'); zlabel('z');
    end
    
    hold off
     light
    
    axis equal
%     view(-270,0);
%     view(-55, 15);
% view([-65 6.2]);
% view([-186 5]);
% view([0 90]);
% view([-25 20]);
view([-45.6 17]);
%     xlim([-2 4.5]); ylim([-5 5]); zlim([-5 5]);
%      xlim([-1.5 1.5]); ylim([-2 2]); zlim([-2 2]);
%        xlim([-1.5 1.5]); ylim([-0.75 0.75]); zlim([-0.75 0.75]);
%         xlim([-1.5 3]); ylim([-3.5 3.5]); zlim([-3.5 3.5]);
%          xlim([-2.2 1]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
%          xlim([-5.5 2.5]); ylim([-2 2]); zlim([-2 2]);
           xlim([-4 5.5]); ylim([-4.5 4.5]); zlim([-4.5 4.5]);
         
         
%     title({"t = " + num2str(t) , "# offenders = " + num2str(length(insides)), "relative distance of worst offender = " + num2str(- min(0,min(rel_dists)))})
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

toc