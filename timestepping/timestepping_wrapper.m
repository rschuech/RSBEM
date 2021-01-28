function [sol, stored_output, trespassers ] = timestepping_wrapper(sol, stored_output, tfinal, Mesh0, Network0, mesh_node_parameters, index_mapping, matrix_props, assembly_input,timestepping_tol, solver)

stored_output_growth_increment = 100; % grow by this amount when needed to avoid slowdown due to growing every timestep

clear textprogressbar output_fun odefun
options = odeset('Refine',1,'RelTol',timestepping_tol.reltol,'AbsTol',timestepping_tol.abstol,'OutputFcn',@output_fun,'Events',@trespassing_events); % be sure to use Refine = 1 or else the output of the derivatives won't work right as-is
y0 = [ Network0.l; reshape(Network0.nodes',Network0.n_nodes*3,1) ; 0*Mesh0(1).refpoints(:,1); zeros(4,1)];

store_network_velocities = false;  % to store velocity of every network node or not (kinematic variables of swimmer are always stored)
n_kinematic_variables = matrix_props.n_cols - matrix_props.n_collocation*3;
start_ind = length(y0) - n_kinematic_variables + 1; % ind of first kinematic variable component i.e. U_1
end_ind = length(y0);

current_time_ind = [];  current_repulsion_ind = [];
% stored_t_y = [];
% is_inside = [];

if isempty(sol)
    tspan = [0 tfinal];
    current_time_ind = 1;  %index for every output timestep
    current_repulsion_ind = 1; % index for every timestep in which there were any network nodes close to or inside swimmer
    
    if store_network_velocities
        stored_output.derivatives = NaN(stored_output_growth_increment, length(y0));
    else
        stored_output.derivatives = NaN(stored_output_growth_increment, n_kinematic_variables);
    end
    
    stored_output.time = NaN(stored_output_growth_increment,1);
    
    stored_output.link_breakage_time = NaN(Network0.n_links,1); % for each link, will be either NaN if link never broke or the time that it broke (and remained broken thereafter)
    
    
    
    stored_output.PE = NaN(stored_output_growth_increment,1);
    
    stored_output.repulsion.time = NaN(stored_output_growth_increment,1);
    stored_output.repulsion.any_flagged = false(stored_output_growth_increment,1);
    stored_output.repulsion.any_inside = false(stored_output_growth_increment,1);
    stored_output.repulsion.network_node_indices = cell(stored_output_growth_increment,1);
    stored_output.repulsion.submesh_index = cell(stored_output_growth_increment,1);
    stored_output.repulsion.distance2mesh = cell(stored_output_growth_increment,1);
    stored_output.repulsion.is_inside = cell(stored_output_growth_increment,1);
    stored_output.repulsion.total_force_mag = zeros(stored_output_growth_increment,1);
    stored_output.repulsion.total_torque_mag = zeros(stored_output_growth_increment,1);
    
    [sol] = solver(@odefun, tspan, y0,options);
    
else
    current_time_ind = length(stored_output.time) ;
    % got rid of + 1 since odextend always seems to repeat the last time
    % value from before and we can assume it will be the same unless
    % there's a discontinuity in the problem
    if stored_output.repulsion.any_flagged(current_time_ind)
        current_repulsion_ind = length(stored_output.repulsion.network_node_indices); %should redo this one and so should overwrite last value with same thing
    else
        current_repulsion_ind = length(stored_output.repulsion.network_node_indices) + 1; % the last entry won't get redone and we should move to next
    end
    [sol] = odextend(sol, @odefun, tfinal,[],options );
    
end

% remove placeholders created during incremental preallocation
ind = find(isnan(stored_output.time),1,'first');
stored_output.time(ind:end) = [];
stored_output.derivatives(ind:end, :) = [];
stored_output.PE(ind:end) = [];
stored_output.repulsion.any_flagged(ind:end) = [];
stored_output.repulsion.any_inside(ind:end) = [];
stored_output.repulsion.total_force_mag(ind:end) = [];
stored_output.repulsion.total_torque_mag(ind:end) = [];

ind = find( cellfun(@isempty, stored_output.repulsion.network_node_indices),1,'first');
stored_output.repulsion.time(ind:end) = [];
stored_output.repulsion.network_node_indices(ind:end) = [];
stored_output.repulsion.submesh_index(ind:end) = [];
stored_output.repulsion.distance2mesh(ind:end) = [];
stored_output.repulsion.is_inside(ind:end) = [];

Mesh = []; Network = [];  trespassers = []; % to make scope across all subfunctions

    function [dydt] = odefun(t,y)
        Network0.E( ~isnan(stored_output.link_breakage_time) ) = 0;  % permanently set E = 0 each time a link breaks
        
        [dydt, Mesh, Network, Repulsion] = derivatives(t,y, Mesh0, Network0, mesh_node_parameters, index_mapping, matrix_props, assembly_input);
        % Network.E == 0 corresponds to links that broke at any time including
        % just now, we will figure out links that just broke by comparing
        % this with stored_output.link_breakage_time that has not been
        % updated yet
        
        if current_time_ind > length(stored_output.time)
            stored_output.time(end+1:end+stored_output_growth_increment,1) = NaN;
            stored_output.derivatives(end+1:end+stored_output_growth_increment, :) = NaN;
            
            stored_output.PE(end+1:end+stored_output_growth_increment,1) = NaN;
            
            stored_output.repulsion.any_flagged(end+1:end+stored_output_growth_increment,1) = false;
            stored_output.repulsion.any_inside(end+1:end+stored_output_growth_increment,1) = false;
            stored_output.repulsion.total_force_mag(end+1:end+stored_output_growth_increment,1) = NaN;
            stored_output.repulsion.total_torque_mag(end+1:end+stored_output_growth_increment,1) = NaN;
            
        end
        
        if current_repulsion_ind > length(stored_output.repulsion.network_node_indices)
            
            stored_output.repulsion.time(end+1:end+stored_output_growth_increment,1) = NaN;
            %below works since entries before
            %end+stored_output_growth_increment automatically filled
            %with [] also
            stored_output.repulsion.network_node_indices{end+stored_output_growth_increment,1} = [];
            stored_output.repulsion.submesh_index{end+stored_output_growth_increment,1} = [];
            stored_output.repulsion.distance2mesh{end+stored_output_growth_increment,1} = [];
            stored_output.repulsion.is_inside{end+stored_output_growth_increment,1} = [];
        end
        
        
        stored_output.time(current_time_ind) = t;
        if store_network_velocities
            stored_output.derivatives(current_time_ind, :) = dydt; % this is shared with timestepping_wrapper and thus also output_fcn
        else
            stored_output.derivatives(current_time_ind, :) = dydt(start_ind : end_ind);
        end
        
        stored_output.link_breakage_time( (Network.E == 0) & isnan(stored_output.link_breakage_time) ) = t; % only store current time value if a link just broke this timestep, i.e. only overwrite NaNs
        
        % stored_t_y = [t; y(end-6:end)];
        %         is_inside = Repulsion.is_inside;
        
        stored_output.PE(current_time_ind,1) = Network.PE;
        stored_output.repulsion.any_flagged(current_time_ind,1) = ~isempty(Repulsion.network_node_indices);
        stored_output.repulsion.any_inside(current_time_ind,1) = any(Repulsion.is_inside);
        
        stored_output.repulsion.total_force_mag(current_time_ind,1) = sqrt(sum( (sum(Repulsion.F(Repulsion.network_node_indices,:), 1)).^2 ))  ;
        [refpoint] = get_rotational_references(Mesh, assembly_input);
        stored_output.repulsion.total_torque_mag(current_time_ind,1) = sqrt(sum(  sum( cross( (Repulsion.x  - refpoint') , Repulsion.F(Repulsion.network_node_indices,:) ), 1).^2  ));
        
        if stored_output.repulsion.any_flagged(current_time_ind,1)
            stored_output.repulsion.time(current_repulsion_ind,1) = t;
            stored_output.repulsion.network_node_indices{current_repulsion_ind,1} = Repulsion.network_node_indices;
            stored_output.repulsion.submesh_index{current_repulsion_ind,1} = Repulsion.mesh_index;
            stored_output.repulsion.is_inside{current_repulsion_ind,1} = Repulsion.is_inside;
            stored_output.repulsion.distance2mesh{current_repulsion_ind,1} = Repulsion.distance2mesh;
        else
            stored_output.repulsion.time(current_repulsion_ind,1) = NaN;
            stored_output.repulsion.network_node_indices{current_repulsion_ind,1} = [];
            stored_output.repulsion.submesh_index{current_repulsion_ind,1} = [];
            stored_output.repulsion.is_inside{current_repulsion_ind,1} = [];
            stored_output.repulsion.distance2mesh{current_repulsion_ind,1} = [];
        end
        
    end

    function status = output_fun(t, y, flag)
        status = 0;
        
        %             current_time_ind = 1;
        %         if isempty(flag)
        
        if stored_output.repulsion.any_inside(current_time_ind)
            %             warning('network node(s) are inside swimmer');
        end
        
        
        if stored_output.repulsion.any_flagged(current_time_ind,1)
            current_repulsion_ind = current_repulsion_ind + 1;
        end
        current_time_ind = current_time_ind + 1;
        
        if find(~cellfun(@isempty,stored_output.repulsion.network_node_indices),1,'last') > sum(stored_output.repulsion.any_flagged)
            warning('wtf')
            pause
        end
        
        %             [stored_t_y   [t; y(end-6:end)]  ]
        %
        % this is a regular successful step, not initialization nor cleanup
        
        %         end
        
        odetpbar(t,y,flag);
        
        
    end


    function [value,isterminal,direction] = trespassing_events(t,y)
        
        % Network_nodes = reshape( y(Network0.n_links + 1 : Network0.n_links + Network0.n_nodes*3) , 3,[])';
        
        % first, filter out any passive network nodes with no links - don't bother worrying about these getting inside, don't bother calculating distances for them
        active = true(Network.n_nodes,1); % active if at least one working link, passive if no active links (all E == 0)
        for i = 1:Network.n_nodes
            active(i) = any( Network.E(  Network.link_members{i}(1,:) ) ~= 0 ); % node is active if it has at least one link with nonzero stiffness
        end
        
%         active_network_nodes = Network.nodes(active,:);
        % any passive nodes default (and remain) as NaN / "not inside" (even though they actually may be inside)
        distance2mesh = NaN(Network.n_nodes,1);
        is_inside = false(Network.n_nodes,1);
        x = NaN(Network.n_nodes,3);
        
        [distance2mesh(active), x(active), xi_eta, mesh_index, element_index, is_inside(active)] = dist2mesh(Network.nodes(active,:), Mesh, index_mapping, assembly_input.repulsion.mindist2);
        % toc
        
        
        inside_threshold = 0.006; % body radius / 100
        % if inside further than this, stop simulation and put node back on surface
        
        % the last thing maps is_inside (0 or 1) to 1 or -1, so we make the distance signed, negative inside
        % then add that to the threshold so that value will be postive if outside or slightly inside, 0 at the threshold, and negative if further inside
        value = inside_threshold + distance2mesh.*(1 - 2*is_inside);
        value(isnan(value)) = 1;  % hopefully this works - if far away outside, or if too far inside - hopefully that will never happen...
        
        isterminal = true(size(value));
        
        direction = [];
        
        trespassers.distance2mesh = distance2mesh;
        trespassers.x = x;
        trespassers.is_inside = is_inside;
        
        
        % [value,isterminal,direction] = myEventsFcn(t,y)
        
        % The output arguments value, isterminal, and direction are vectors whose ith element corresponds to the ith event:
        %
        % value(i) is a mathematical expression describing the ith event. An event occurs when value(i) is equal to zero.
        %
        % isterminal(i) = 1 if the integration is to terminate when the ith event occurs. Otherwise, it is 0.
        %
        % direction(i) = 0 if all zeros are to be located (the default). A value of +1 locates only zeros where the event function is increasing,
        %     and -1 locates only zeros where the event function is decreasing. Specify direction = [] to use the default value of 0 for all events.
        
        
    end





end
