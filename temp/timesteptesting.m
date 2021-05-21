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
    if isfield(sol,'xe') && ~isempty(sol.xe) && sol.xe(end) > last_trespass_time % solver stopped due to an event, not reaching tfinal
        last_trespass_time = sol.xe(end);
        disp('trespassing event')
        sol.xe(end)
        
        %       Network_temp = Network0;
        %        Network_temp.E( ~isnan(stored_output.link_breakage_time) ) = 0;  % permanently set E = 0 each time a link breaks
        % the hell with it, the links don't matter for just figuring out where the
        % network nodes are relative to the swimmer when the event occurs
        [Mesh, Network, mesh_node_parameters, y_swimmer, r] = update_Mesh_Network(Mesh0, Network0, sol.xe(end), sol.ye(:,end), mesh_node_parameters,index_mapping, assembly_input);
        
        [distance2mesh, x, xi_eta, mesh_index, element_index, is_inside] = dist2mesh(Network.nodes, Mesh, index_mapping, assembly_input.repulsion.mindist2);
        
        figure(345);
        [s,e] = plot_mesh(Mesh(1),2); s.FaceAlpha = 0.5; hold on
        pl = plot3(Network.nodes(is_inside,1),Network.nodes(is_inside,2),Network.nodes(is_inside,3),'ko','MarkerFaceColor','k');
        
        
        
        
        %         Network.l = y(1:Network.n_links);
        network_nodes = reshape( sol.ye(Network.n_links + 1 : Network.n_links + Network.n_nodes*3, end) , 3,[])';
        network_nodes( is_inside,: ) = x(is_inside, :);  % put them back on the closest points on the mesh
        % here we are putting all far-inside nodes back on the mesh, even nodes
        % that don't have links - we don't really care about those, but why not fix
        % them too
        delete(pl);
         pl = plot3(network_nodes(is_inside,1),network_nodes(is_inside,2),network_nodes(is_inside,3),'ko','MarkerFaceColor','k');
         
        network_nodes =  reshape(network_nodes', Network.n_nodes*3,1);
        
        y0 = sol.ye(:,end);
        
        y0( Network.n_links + 1 : Network.n_links + Network.n_nodes*3  ) = network_nodes;
        
    else % solver finished normally at the final time
        c = c + 1;  % go to next tfinal
        y0 = sol.y(:,end);
    end
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
    
    
    % end
    %     save('c:\Users\rudi\Desktop\RD\temp\broken.mat','sol');
end