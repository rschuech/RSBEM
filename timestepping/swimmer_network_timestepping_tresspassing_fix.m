%%
Mesh0 = Mesh;  Network0 = Network;
 rng(0);
Network0.nodes = Network0.nodes + rand(size(Network0.nodes)) * 1E-16;
% Network0.E(:) = 0;
 %%
y0_0 = [ Network0.l; reshape(Network0.nodes',Network0.n_nodes*3,1) ; 0*Mesh0(1).refpoints(:,1); zeros(4,1)];
assembly_input.performance.verbose = false;
assembly_input.Tail.motor_torque = @(t)motor_torque(t,torque0,1000);
clear textprogressbar
clear matrix_assembly_mex
solver = @ode45;

sol = []; stored_output = [];
last_trespass_time = -Inf;
y0 = y0_0;
% 'RelTol',5E-4,'AbsTol',5E-4,  originally
% defaults are 1E-3, 1E-6
% timestepping_tol.reltol = 1E-6;
% timestepping_tol.abstol = 1E-10;
% timestepping_tol.reltol = 5E-4;
% timestepping_tol.abstol = 5E-4;
timestepping_tol.reltol = 1E-3 /10 * Inf  ;
timestepping_tol.abstol = 1E-8  /10 * Inf ;
% timestepping_tol.reltol = 5E-4;
% timestepping_tol.abstol = 5E-4;
timestepping_tol.InitialStep = [];
timestepping_tol.MaxStep = 1E-3 ;
% [sol, stored_output, trespassers ] = timestepping_wrapper([], [],y0, 0.01, Mesh0, Network0, mesh_node_parameters, index_mapping, matrix_props, assembly_input,timestepping_tol, solver);

clear textprogressbar
% tfinals = 0.005:0.005:2;
tfinals = 0.02:0.01:1.4;
c = 1;
n_events_prior = 0;
total_time = tic;
%%
while isempty(sol) || sol.x(end) < tfinals(end)
    
  
    
    if ~isempty(sol) && sol.x(end) >= tfinals(c)
        c = c + 1;
        continue
    end
    
    [sol, stored_output , trespassers] = timestepping_wrapper(sol, stored_output,y0, tfinals(c), Mesh0, Network0, mesh_node_parameters, index_mapping, matrix_props, assembly_input,timestepping_tol, solver);
    % if any(trespassers.is_inside)
    %   if any(trespassers.distance2mesh(trespassers.is_inside) >= 0.006 )
    %       pause
%     if isfield(sol,'xe') && ~isempty(sol.xe) && sol.xe(end) > last_trespass_time % solver stopped due to an event, not reaching tfinal
%         last_trespass_time = sol.xe(end);
%         disp('trespassing event')
%         sol.xe(end)
         if isfield(sol,'xe') && ~isempty(sol.xe) && length(sol.xe) > n_events_prior
             1;
             n_events_prior = length(sol.xe);
         end
        %       Network_temp = Network0;
        %        Network_temp.E( ~isnan(stored_output.link_breakage_time) ) = 0;  % permanently set E = 0 each time a link breaks
        % the hell with it, the links don't matter for just figuring out where the
        % network nodes are relative to the swimmer when the event occurs
      
        % this is where we put nodes back for the hell of it
%         [Mesh, Network, mesh_node_parameters, y_swimmer, r] = update_Mesh_Network(Mesh0, Network0, sol.x(end), sol.y(:,end), mesh_node_parameters,index_mapping, assembly_input);
%         [distance2mesh, x, xi_eta, mesh_index, element_index, is_inside] = dist2mesh(Network.nodes, Mesh, index_mapping, assembly_input.repulsion.mindist2);
%         
%         figure(345);  clf;
%         [s,e] = plot_mesh(Mesh(1),2); s.FaceAlpha = 0.5; hold on
%         pl = plot3(Network.nodes(is_inside,1),Network.nodes(is_inside,2),Network.nodes(is_inside,3),'ko','MarkerFaceColor','k');
%         
%         network_nodes = reshape( sol.y(Network.n_links + 1 : Network.n_links + Network.n_nodes*3, end) , 3,[])';
%         network_nodes( is_inside,: ) = x(is_inside, :);  % put them back on the closest points on the mesh
%         % here we are putting all far-inside nodes back on the mesh, even nodes
%         % that don't have links - we don't really care about those, but why not fix
%         % them too
%         delete(pl);
%          pl = plot3(network_nodes(is_inside,1),network_nodes(is_inside,2),network_nodes(is_inside,3),'ko','MarkerFaceColor','k');
%          
%         network_nodes =  reshape(network_nodes', Network.n_nodes*3,1);
        
        y0 = sol.y(:,end);
        
%         y0( Network.n_links + 1 : Network.n_links + Network.n_nodes*3  ) = network_nodes;
        
%     else % solver finished normally at the final time
%         c = c + 1;  % go to next tfinal
%         y0 = sol.y(:,end);
% end
 


    %         [sol, stored_output , trespassers] = timestepping_wrapper([sol], [stored_output],y0_new, tfinal, Mesh0, Network0, mesh_node_parameters, index_mapping, matrix_props, assembly_input,timestepping_tol, solver);
    
    
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
    
    if sol.x(end) == 0.05
%          save('c:\Users\rudi\Desktop\RD\temp\up_to_0.05.mat','sol','stored_output');
    end
    % end
    %     save('c:\Users\rudi\Desktop\RD\temp\broken.mat','sol');
      "link breakage body 0 tail 0, E 10 eta 500, maxstep 1E-3, iterative refinement, perturbation"
      sol.x(end)
      toc(total_time) / 3600
      
end


return

sqrt(assembly_input.accuracy.mesh.eps2) % 0.4E-4
input.accuracy.mesh.integration_tol.stokeslet.abstol % 0.2
sol.extdata.options