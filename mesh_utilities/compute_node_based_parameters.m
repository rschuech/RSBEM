function [node_parameters, index_mapping, Mesh] = compute_node_based_parameters(Mesh,index_mapping, parent_topologies, coincident_submeshes, BC_type)
% creates node_parameters containing node-specific parameters e.g. avg normal, avg tangents, solid angle
% also add to index_mappings: global2local conversions as well as mappings for unknown u at free slip nodes
% also add to Mesh:  solid angle (ordinarily we'd only need an overall solid angle for each global node but if one decides to only eliminate DL on some but
% not all entities that share a node, we need solid angles for each isolated entity considered alone

%#codegen

%% go from a global node index to all of its instances as a local node in an element of a submesh
% while we're here, also initialize solid angle fields of Mesh

% n_global_nodes = length(index_mapping.global_node2local_node);
n_global_nodes = max( vertcat( index_mapping.local_node2global_node{:} ) );

index_mapping.global_node2local_node = cell( n_global_nodes , 1 );
% index_mapping.global_u2global_node = NaN(  % difficult to initialize since we don't know exactly how many unknown u there are yet I guess
    
for submesh_ind = 1:length(index_mapping.local_node2global_node)  %submeshes
    
    for local_node_ind = 1:Mesh(submesh_ind).n_nodes
        global_ind = index_mapping.local_node2global_node{submesh_ind}(local_node_ind);
        [element_inds, position_inds] = find( Mesh(submesh_ind).elements == local_node_ind );
        
        % global2local is a n_global_nodes cell array, each cell is a matrix with rows = submesh, local node ind, element ind, position in element
        index_mapping.global_node2local_node{global_ind}(end+1:end+size(element_inds,1),:) = [repmat([submesh_ind, local_node_ind],size(element_inds,1),1), element_inds, position_inds];
    end
    
    Mesh(submesh_ind).solid_angle = NaN(Mesh(submesh_ind).n_nodes,1);
end




node_parameters.node_type = NaN(n_global_nodes,1); % 1 for vertex node, 2 for edge (midpoint) node on the curved triangular elements
% node_parameters.solid_angle = NaN(n_global_nodes,1);
node_parameters.normals_avg = NaN(n_global_nodes,3);
node_parameters.tangents_avg = NaN(n_global_nodes,3,2);
node_parameters.BC_type = NaN(n_global_nodes, 1);
index_mapping.global_node2global_u = zeros(n_global_nodes,1);

names = [Mesh.name];
topologies = strings(length(names),1); BC_types = NaN(length(names),1);
for n = 1:length(names)
    topologies(n) = parent_topologies.(names(n));
    BC_types(n) = BC_type.(names(n));
end

closed_submeshes = find(topologies == "closed"); % submesh inds of closed submeshes, since we don't want to account for open submesh parents for solid angle


entities = [];
for i = 1:length(coincident_submeshes)
    for j = 1:length(coincident_submeshes{i})
        entities{end+1} = [coincident_submeshes{i}{j}];
    end
end

% last_unknown_u_global_ind = 0; % will increment later as we add global nodes with unknown u due to free slip BC

for global_node_ind = 1:n_global_nodes
    
    
    n_local_entries = size(index_mapping.global_node2local_node{global_node_ind},1);
    normals = NaN(n_local_entries,3);
    angles = NaN(n_local_entries,1);
    
    
    % decide what BC to impose for this global node
    if any( BC_types( unique(index_mapping.global_node2local_node{global_node_ind}(:,1)) ) == 1 ) % is the node a member of at least one no-slip submesh?
        
        node_parameters.BC_type(global_node_ind, 1) = 1; % arbitrarily decide that the node is then no-slip (even if also a member of free-slip submeshes)
%         index_mapping.global_node2global_u(global_node_ind) = 0; % this global node doesn't have an unknown u
    else
        node_parameters.BC_type(global_node_ind, 1) = 2; % if node is only a member of free-slip submeshes, it must obviously be free-slip
        index_mapping.global_node2global_u(global_node_ind) = max( index_mapping.global_node2global_u ) + 1;
        index_mapping.global_u2global_node( index_mapping.global_node2global_u(global_node_ind) , 1 ) = global_node_ind;
    end
    
    
    for local_entry = 1:n_local_entries % each time this node appears in an element, over all submeshes
        
        
        parent_submesh =  index_mapping.global_node2local_node{global_node_ind}(local_entry,1);
        local_node_ind = index_mapping.global_node2local_node{global_node_ind}(local_entry,2);
        
        
        containing_element = index_mapping.global_node2local_node{global_node_ind}(local_entry,3); % actual local element index
        node_position = index_mapping.global_node2local_node{global_node_ind}(local_entry,4); % node position within element, 1 - 6
        
        shape_parameters = Mesh(parent_submesh).shape_parameters(containing_element,:)'; % alpha beta gamma
        
        node_coords = Mesh(parent_submesh).nodes(Mesh(parent_submesh).elements(containing_element,:),:); % coords of the 6 nodes in this element
        
        switch node_position % see Pozrikidis page 123 Fig 5.3.1
            % rejigger alpha, beta, gamma so that the corner vertex (nodes 1,2,3) is the "first" one, at xi = eta = 0, so that computed tangents are along the
            % sides connected to this vert
            case 1
                xi_eta = [0 0]; % node_coords = node_coords([1 2 3 4 5 6],:); shape_parameters = [shape_parameters(1);  shape_parameters(2);  shape_parameters(3) ];
            case 2
                xi_eta = [0 0];  node_coords = node_coords([2 3 1 5 6 4],:);  shape_parameters = [1 - shape_parameters(3);   shape_parameters(2);   1 - shape_parameters(1) ];
            case 3
                xi_eta = [0 0];  node_coords = node_coords([3 1 2 6 4 5],:);  shape_parameters = [1 - shape_parameters(2);  1 - shape_parameters(1);  shape_parameters(3) ];
                
                % don't need to do any rejiggering for midpoint nodes (nodes 4,5,6)
            case 4
                xi_eta = [shape_parameters(1) 0];
            case 5
                xi_eta = [shape_parameters(3), 1 - shape_parameters(3)];
            case 6
                xi_eta = [0 shape_parameters(2)];
        end
        
        
        % below stores normals based on all elements that each node is a member of
        % not permanently storing since not sure these would ever be useful beyond calculating avg normals
        switch node_position
            case {1 2 3} % triangle verts
                [~, ~, ~,  normals(local_entry,:), tangents_along_edge] = T6interp(node_coords,xi_eta(1),xi_eta(2),shape_parameters);
                node_parameters.node_type(global_node_ind) = 1; % vertex node (will get overwritten with same value multiple times, who cares)
                angles(local_entry,1) = acos(dot( tangents_along_edge(:,1) , tangents_along_edge(:,2 )));
            case {4 5 6} % triangle edge midpoints
                [~, ~, ~,  normals(local_entry,:), ~] = T6interp(node_coords,xi_eta(1),xi_eta(2),shape_parameters);
                node_parameters.node_type(global_node_ind) = 2; % midpoint node (will get overwritten with same value multiple times, who cares)
                angles(local_entry,1) = NaN;
        end
        
        
        
    end % each instance of this global node across submeshes, elements
    
    
    
    if size( normals , 1 ) > 2 % if node belongs to at least 3 elements, take weighted avg of normals, weighted with element angles
        temp = sum( normals .* angles, 1 ) / sum(angles);
    else % node belongs to 1 or 2 elements, simply take average of normals (avg normal = single normal for 1 element case)
        % note:  only way a node can belong to 2 elements on a closed surface is a midpoint node, can't belong to just 1 element; only way a node can belong to 2 elements on a sheet is midpoint node
        % or a node on the edge of the sheet; only way a node can belong to 1 element on a sheet is for a corner of the sheet
        % I *think* even if there is a cusp along the edge of a sheet (pretty bizarre, we're saying the cusp is within an element vs logically occuring along
        % the edges of adjacent elements), the avg normal should be a simple unweighted average, the same as a midpoint along a ridge in a closed surface mesh
        % here we are considering elements on each submesh separately, so a node shared by multiple submeshes will have multiple avg normals etc
        temp = mean( normals , 1 ) ;
    end
    node_parameters.normals_avg(global_node_ind,:) = temp / sqrt(sum(temp.^2 , 2));
    
    
    % there's nothing "wrong" with along-element-edge tangents that can be output from T6interp but if we want an "avg" tangent at a corner point, it
    % isn't straightforward to calculate it from those tangents (e.g. for very symmetric element arrangements, an average of all tangents can yield
    % nearly zero due to cancellations) so instead, compute tangents from just the avg normal
    tangents = eye(3) - node_parameters.normals_avg(global_node_ind,:)' * node_parameters.normals_avg(global_node_ind,:); % Dave, or common knowledge?...
    tangents = tangents ./ sqrt(sum(tangents.^2,2));
    % the above yields 3 tangent vectors.  choose the pair that are most orthogonal to each other for giggles (as long as the angle between them isn't
    % nearly zero, it's fine for the purposes of imposing a free-slip BC)
    dots = [dot(tangents(3,:),tangents(2,:)); dot(tangents(1,:),tangents(3,:)); dot(tangents(2,:),tangents(1,:))];
    [~,ind] = min(abs(dots));
    tangents(ind,:) = []; % delete the extra tangent vector
    node_parameters.tangents_avg(global_node_ind,:,:) = tangents'; % 3rd dim is for the two tangent vectors
    
    % solid angle is calculated according to general method in following, p. 24
    % Topping, B. H. V., J. Muylle, R. Putanowicz, and P. Ivanyi. Finite Element Mesh Generation. Paul & Company Pub Consortium, 2004.
    
    % phi is the angle between normals of adjacent elements (for our curved meshes, at an edge or corner, the value of the calculated normal at the
    % discontinuity will depend on which element it is assumed to be a part of)
    
    
    
    %     inputs.coincident_submeshes = {"Body" , "Tail"  };
    
     % global2local is a n_global_nodes cell array, each cell is a matrix with rows = submesh, local node ind, element ind, position in element
    filter = ismember(index_mapping.global_node2local_node{global_node_ind}(:,1), closed_submeshes); % global node instances that are on closed submeshes
    filtered_global2local = index_mapping.global_node2local_node{global_node_ind}(filter,:);
    
    
%     entities = [coincident_submeshes{:}]; % this doesn't seem to work, can't figure out a one liner, so doing it the n00b way above
    for e = 1:length(entities)
        [~,submesh_inds] = ismember( entities{e}  , [Mesh.name]  );
        is_open = false(length(entities{e}),1);
        
        for s = 1:length(entities{e})
            if parent_topologies.(entities{e}(s)) == "open"  % if any member of a parent claims its parent is open, we can call the entire entity open
                is_open(s) = true;
                % this global node is only part of open submeshes, which by definition take up (nearly) zero fluid volume, so phi = 4*pi (?)
                %           submesh_ind = find([Mesh.name] == entities{e}(s));
                local_node_ind = index_mapping.global_node2local_node{global_node_ind}(find(index_mapping.global_node2local_node{global_node_ind}(:,1) == submesh_inds(s) , 1),2);
                % shouldn't matter which instance of this node we find in this submesh, so just find first one
                Mesh(submesh_ind(s)).solid_angle(local_node_ind) = 4*pi; % current thinking is for an infinitely thin (or eps*2 width) sheet, solid angle should be a full sphere,
                %                 no fluid volume obscured by geometry?  Doesn't matter in the end since solid angle for sheets cancels out in the BIE.
    
            end
          
        end
        
        if ~ismember(sum(is_open) , [length(is_open), 0]) % all submeshes of a geometric entity should claim to either be part of a closed or open surface
            error("all submeshes of a geometric entity should claim to either be part of a closed or open surface - check inputs structure");
        end
        
        if all(is_open)
            continue
        end
        
        %            [~,inds] = ismember( entities{e}  , [Mesh.name]  );
        filter = ismember(index_mapping.global_node2local_node{global_node_ind}(:,1), submesh_inds); % global node instances that are on this closed entity
        filtered_global2local = index_mapping.global_node2local_node{global_node_ind}(filter,:); % only contains instances of current global node on submeshes
        % that are part of current closed geometric entity
        
        phi = [];
        
        
        containing_elements = filtered_global2local(:,3); % actual local element index
        parent_submeshes = filtered_global2local(:,1);
        
        for v1 = 1:length(containing_elements) - 1 % elements containing node of interest, v1 goes with "first" potential neighbor under interrogation
            % global indices of the 6 nodes belonging to element v1:
            global_nodes1 = index_mapping.local_node2global_node{parent_submeshes(v1)}( Mesh(parent_submeshes(v1)).elements(containing_elements(v1),:) );
%             global_nodes1(filtered_global2local(v1,4)) = NaN; % don't want to consider the current global node when looking for adjacent elements
            for v2 = v1+1:length(containing_elements) % v2 goes with "second" potential neighbor under interrogation
                global_nodes2 = index_mapping.local_node2global_node{parent_submeshes(v2)}( Mesh(parent_submeshes(v2)).elements(containing_elements(v2),:) );
%                 global_nodes2(filtered_global2local(v1,4)) = NaN; % don't want to consider the current global node when looking for adjacent elements
                % NaNing out doesn't seem to make sense, since current global node must be one of the shared nodes between adjacent elements, if they are
                % adjacent
                %                 if ~isempty( intersect( global_nodes1(1:3) , global_nodes2(1:3) ) ) % elements v1 and v2 are adjacent neighbors sharing an edge
                %                     neighbors = [containing_elements(v1) containing_elements(v2)]; %local element indices of neighbor 1 and neighbor 2
                %                     Mesh(i).elem_neighbors{vi}(end+1,:) = neighbors;
                [~,i1,i2] = intersect(  global_nodes1(1:3)  ,  global_nodes2(1:3) );  % global inds of the 2 shared corner vertices of adjacent elements, if they are adjacent
                if length(i1) == 2 % these two elements are in fact adjacent
                    % i1, i2 tell us the positions of the shared verts inside the two elements.  first get the position ind (1,2,or 3) that is NOT
                    % shared in each.  then get the local node inds of these verts.  then get their coordinates
                    vert1 =  Mesh(parent_submeshes(v1)).nodes( Mesh(parent_submeshes(v1)).elements(containing_elements(v1), setdiff(1:3,i1)  ) , : );
                    vert2 =  Mesh(parent_submeshes(v2)).nodes( Mesh(parent_submeshes(v2)).elements(containing_elements(v2), setdiff(1:3,i2)  ) , : );
                    
                    %                         https://math.stackexchange.com/questions/1790118/determine-the-concavity-of-an-edge
                    if dot( (vert2 - vert1) , normals(v1,:) ) <= 0
                        phi_sign = -1;
                    else
                        phi_sign = 1;
                    end
                    phi(end+1,1) = phi_sign * acos( min(max(dot(normals( v1,:)  ,  normals( v2,:) ) , -1),1)  );
                    % occasionally roundoff error yields a dot product outside [-1, 1], so limit dot prod to this range to avoid complex acos
                    %
                end
            end
        end
        
        
        
        %     if ~isempty(closed_submeshes) % this global node is part of at least one closed body
        for s = submesh_inds
            local_node_ind = index_mapping.global_node2local_node{global_node_ind}(find(index_mapping.global_node2local_node{global_node_ind}(:,1) == s , 1),2);
            % shouldn't matter which instance of this node we find in this submesh, so just find first one
            switch length(phi)
                case 1 % this node is a member of exactly 2 neighboring elements, must be a midpoint node
                    Mesh(s).solid_angle(local_node_ind) = 2*pi - phi*2; % not positive about this but seems right for e.g. a midpoint along a ridge
                    % (as unlikely as a discontinuity at a midpoint node would seem), so hopefully generally correct
                otherwise % this node is a member of at least 3 neighboring elements, must be a corner vert on a closed surface
                    %                     Topping 2015 section 2.4.5
                    Mesh(s).solid_angle(local_node_ind,1) = 2*pi - sum(phi);
            end
        end
        %     else % this global node is only part of open submeshes, which by definition take up (nearly) zero fluid volume, so phi = 4*pi (?)
        %         node_parameters.solid_angle(global_node_ind,1) = 4*pi; % current thinking is for an infinitely thin (or eps*2 width) sheet, solid angle should be a full sphere,
        %         %                 no fluid volume obscured by geometry?  Doesn't matter in the end since solid angle for sheets cancels out in the BIE.
        %     end
        
    end % closed entities made of possibly multiple submeshes
    
end % global nodes
