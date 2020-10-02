
function [Mesh, node_parameters] = store_mesh_constants_mex(Mesh,index_mapping, input)

% node_parameters stores data pertaining to each node in a global sense, in contrast to the data attached to each submesh in Mesh, that mostly pertains to
% elements (which are never shared across submeshes)


%#codegen


n_global_nodes = length(index_mapping.global2local.node);

node_parameters.node_type = NaN(n_global_nodes,1); % 1 for vertex node, 2 for edge (midpoint) node on the curved triangular elements

for i = 1:length(Mesh)  %loop over submeshes in case there is more than one
    
    %all integrals below are gratuitously adaptive, using tolerances in input struct
    
    centroid = NaN(3,size(Mesh(i).elements,1)); %centroid of each element

    element_parameters = NaN(size(Mesh(i).elements,1),3);


    %     parfor (elem_i = 1:size(Mesh(i).elements,1), input.performance.nthreads)
    for elem_i = 1:size(Mesh(i).elements,1) %parfor broke when normals and tangents calc was added, prolly not worth fixing
        %elem_i / size(elements,1)
        node_coords = Mesh(i).nodes(Mesh(i).elements(elem_i,:),:); % coords of the 6 nodes in this element
        % see Pozrikidis page 123 for what alpha, beta, gamma mean:
        alpha = 1/(1+norm(node_coords(4,:)-node_coords(2,:))/norm(node_coords(4,:)-node_coords(1,:))); %corresponds to location of edge node # 4
        beta = 1/(1+norm(node_coords(6,:)-node_coords(3,:))/norm(node_coords(6,:)-node_coords(1,:))); % corresponds to location of edge node # 6
        gamma = 1/(1+norm(node_coords(5,:)-node_coords(2,:))/norm(node_coords(5,:)-node_coords(3,:))); % corresponds to location of edge node # 5
        
        element_parameters(elem_i,:) = [alpha beta gamma];
        
        centroid(:,elem_i) = T6interp(node_coords,1/3,1/3,element_parameters(elem_i,:)); %coords of element centroid
   
    end % elements
    
    Mesh(i).element_parameters = element_parameters;
    Mesh(i).centroid = centroid;
  
    Mesh(i).area = ones(Mesh(i).n_elements,1); %need to temporarily initialize this as something since abstol for integration is always scaled by element area, but we're trying to calculate element area....
    %'starting area'
    [VL_unity] = integrate_unity(Mesh(i), input,'area'); %compute integrals of hS * phi
    %    'area done'
    Mesh(i).area = sum(VL_unity,2);  %area is just integral of 1 * hS * phi so can just add up all the hS * phi integrals for each element
    %     'starting flux'
    VL_flux = integrate_flux(Mesh(i), input, 'volume');  %compute integrals of F dot n * hS * phi ala divergence theorem, F hardcoded as [1/3 x  1/3 y  1/3 z] but other F can be added later
    %     'flux done'
    Mesh(i).Volume = sum(VL_flux(:));  %not sure integrals for each individual element are very useful here, but the whole thing is the total volume
    
    %     'starting centroid'
    [VL_centroid] = centroid_integrals(Mesh(i), input);  %output:  row is elem, cols are phi = 1..6 for x, y, z centroid integrals
    %    'centroid done'
    Mesh(i).Centroid(1,1) = sum(sum(VL_centroid(:,1:6))) / Mesh(i).Volume; %not sure if finer grained sums mean anything, but complete sum over all elements is integral (x dV)  and then normalize by total volume
    Mesh(i).Centroid(2,1) = sum(sum(VL_centroid(:,7:12))) / Mesh(i).Volume;
    Mesh(i).Centroid(3,1) = sum(sum(VL_centroid(:,13:18))) / Mesh(i).Volume;
    
end % submeshes



node_parameters.solid_angle = NaN(n_global_nodes,1);

node_parameters.normals_avg = NaN(n_global_nodes,3);
node_parameters.tangents_avg = NaN(n_global_nodes,3,2);


topologies = []; BC_types = [];
names = [Mesh.name];
for n = 1:length(names)
    topologies(n) = input.parent_topology.(names(n));
    BC_types(n) = input.BC_type.(names(n));
end

last_unknown_u_global_ind = 0; % will increment later as we add global nodes with unknown u due to free slip BC

for global_node_ind = 1:n_global_nodes

   
    n_local_entries = size(index_mapping.global2local.node{global_node_ind},1);
    normals = NaN(n_local_entries,3);
    angles = NaN(n_local_entries,1);
    
    
    % decide what BC to impose for this global node
    if any( BC_types( unique(index_mapping.global2local.node{global_node_ind}(:,1)) ) == 1 ) % is the node a member of at least one no-slip submesh?
        
        node_parameters.BC_type(global_node_ind) = 1; % arbitrarily decide that the node is then no-slip (even if also a member of free-slip submeshes)
    else
        node_parameters.BC_type(global_node_ind) = 2; % if node is only a member of free-slip submeshes, it should obviously be free-slip
    end
    
    
    for local_entry = 1:size(index_mapping.global2local.node{global_node_ind},1) % each time this node appears in an element, over all submeshes
        
      
        parent_submesh =  index_mapping.global2local.node{global_node_ind}(local_entry,1);
        local_node_ind = index_mapping.global2local.node{global_node_ind}(local_entry,2);
      
        
        containing_element = index_mapping.global2local.node{global_node_ind}(local_entry,3); % actual local element index
        node_position = index_mapping.global2local.node{global_node_ind}(local_entry,4); % node position within element, 1 - 6

        element_parameters = Mesh(parent_submesh).element_parameters(containing_element,:)'; % alpha beta gamma
        
        node_coords = Mesh(parent_submesh).nodes(Mesh(parent_submesh).elements(containing_element,:),:); % coords of the 6 nodes in this element
        
        switch node_position % see Pozrikidis page 123 Fig 5.3.1
            % rejigger alpha, beta, gamma so that the corner vertex (nodes 1,2,3) is the "first" one, at xi = eta = 0, so that computed tangents are along the
            % sides connected to this vert
            case 1
                xi_eta = [0 0]; % node_coords = node_coords([1 2 3 4 5 6],:); element_parameters = [element_parameters(1);  element_parameters(2);  element_parameters(3) ];
            case 2
                xi_eta = [0 0];  node_coords = node_coords([2 3 1 5 6 4],:);  element_parameters = [1 - element_parameters(3);   element_parameters(2);   1 - element_parameters(1) ];
            case 3
                xi_eta = [0 0];  node_coords = node_coords([3 1 2 6 4 5],:);  element_parameters = [1 - element_parameters(2);  1 - element_parameters(1);  element_parameters(3) ];
                
                % don't need to do any rejiggering for midpoint nodes (nodes 4,5,6)
            case 4
                xi_eta = [element_parameters(1) 0];
            case 5
                xi_eta = [element_parameters(3), 1 - element_parameters(3)];
            case 6
                xi_eta = [0 element_parameters(2)];
        end
        
        
        % below stores normals based on all elements that each node is a member of
        % not permanently storing since not sure these would ever be useful beyond calculating avg normals
        switch node_position
            case {1 2 3} % triangle verts
                [~, ~, ~,  normals(local_entry,:), tangents_along_edge] = T6interp(node_coords,xi_eta(1),xi_eta(2),element_parameters);
                node_parameters.node_type(global_node_ind) = 1; % vertex node (will get overwritten with same value multiple times, who cares) 
                angles(local_entry,1) = acos(dot( tangents_along_edge(:,1) , tangents_along_edge(:,2 )));
            case {4 5 6} % triangle edge midpoints
                [~, ~, ~,  normals(local_entry,:), ~] = T6interp(node_coords,xi_eta(1),xi_eta(2),element_parameters);
                node_parameters.node_type(global_node_ind) = 2; % midpoint node (will get overwritten with same value multiple times, who cares) 
                angles(local_entry,1) = NaN;
        end
        
        
        if  node_parameters.BC_type(global_node_ind) == 1  % no slip
            index_mapping.local2global(parent_submesh).unknown_u( local_node_ind, 1  ) = 0; % global unknown_u index set to zero for no-slip nodes
        else % free slip
            if local_entry == 1 % first time we've found this free slip global node in any elements here
                last_unknown_u_global_ind = last_unknown_u_global_ind + 1;
            end
            index_mapping.local2global(parent_submesh).unknown_u( local_node_ind, 1  ) = last_unknown_u_global_ind;
            index_mapping.global2local.unknown_u{last_unknown_u_global_ind}(end+1,:) = [parent_submesh, local_node_ind];
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
    node_parameters.normal_avg(global_node_ind,:) = temp / sqrt(sum(temp.^2 , 2));
    
    
    % there's nothing "wrong" with along-element-edge tangents that can be output from T6interp but if we want an "avg" tangent at a corner point, it
    % isn't straightforward to calculate it from those tangents (e.g. for very symmetric element arrangements, an average of all tangents can yield
    % nearly zero due to cancellations) so instead, compute tangents from just the avg normal
    tangents = eye(3) - node_parameters.normal_avg(global_node_ind,:)' * node_parameters.normal_avg(global_node_ind,:); % Dave, or common knowledge?...
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
    
    phi = [];
    
    
    closed_submeshes = find(topologies == "closed"); % submesh inds of closed submeshes, since we don't want to account for open submesh parents for solid angle
    
    filter = ismember([index_mapping.global2local.node{global_node_ind}(:,1)], closed_submeshes); % global node instances that are on closed submeshes
    filtered_global2local = index_mapping.global2local.node{global_node_ind}(filter,:);
    
    
    containing_elements = filtered_global2local(:,3); % actual local element index
    parent_submeshes = filtered_global2local(:,1);
    
    for v1 = 1:length(containing_elements) - 1 % elements containing node of interest, v1 goes with "first" potential neighbor under interrogation
        % global indices of the 6 nodes belonging to element v1:
        global_nodes1 = index_mapping.local2global(parent_submeshes(v1)).node( Mesh(parent_submeshes(v1)).elements(containing_elements(v1)) );
        global_nodes1(filtered_global2local(v1,4)) = NaN; % don't want to consider the current global node when looking for adjacent elements
        for v2 = v1+1:length(containing_elements) % v2 goes with "second" potential neighbor under interrogation
            global_nodes2 = index_mapping.local2global(parent_submeshes(v2)).node( Mesh(parent_submeshes(v2)).elements(containing_elements(v2)) );
            global_nodes2(filtered_global2local(v1,4)) = NaN; % don't want to consider the current global node when looking for adjacent elements
            
            %                 if ~isempty( intersect( global_nodes1(1:3) , global_nodes2(1:3) ) ) % elements v1 and v2 are adjacent neighbors sharing an edge
            %                     neighbors = [containing_elements(v1) containing_elements(v2)]; %local element indices of neighbor 1 and neighbor 2
            %                     Mesh(i).elem_neighbors{vi}(end+1,:) = neighbors;
            [~,i1,i2] = intersect(  global_nodes1(1:3)  ,  global_nodes2(1:3) );  % global inds of the 2 shared corner vertices of adjacent elements, if they are adjacent
            if ~isempty(i1) % these two elements are in fact adjacent
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
    
    
    
    if ~isempty(closed_submeshes) % this global node is part of at least one closed body
        
        switch length(phi)
            case 1 % this node is a member of exactly 2 neighboring elements, must be a midpoint node
                node_parameters.solid_angle(global_node_ind,1) = 2*pi - phi*2; % not positive about this but seems right for e.g. a midpoint along a ridge
                % (as unlikely as a discontinuity at a midpoint node would seem), so hopefully generally correct
            otherwise % this node is a member of at least 3 neighboring elements, must be a corner vert on a closed surface
                %                     Topping 2015 section 2.4.5
                node_parameters.solid_angle(global_node_ind,1) = 2*pi - sum(phi);
        end
        
    else % this global node is only part of open submeshes, which by definition take up (nearly) zero fluid volume, so phi = 4*pi (?)
        node_parameters.solid_angle(global_node_ind,1) = 4*pi; % current thinking is for an infinitely thin (or eps*2 width) sheet, solid angle should be a full sphere,
        %                 no fluid volume obscured by geometry?  Doesn't matter in the end since solid angle for sheets cancels out in the BIE.
    end
    
    
    
end % global nodes
