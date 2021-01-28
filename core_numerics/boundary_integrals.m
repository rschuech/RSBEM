function [BI_coll, RHS_coll, BI_coll0, RHS_coll0] = boundary_integrals(x0_parameters, global_node_ind, Mesh, mesh_node_parameters, Network, Repulsion, matrix_props, index_mapping, BI_parameters , firstrun)

% for field vel calc, input some placeholder global_node_ind like NaN since the field point (e.g. network node) doesn't have a global mesh node ind, it is
% higher than the bounds of the global mesh node inds

% perhaps input actual global node ind, check if within or above the bounds for the meshes, if above, then do shat for network node

% need to add loop over network nodes at end of this

% this loops over collocation points
%         for global_node_ind = 1:matrix_props.n_collocation    %randperm(matrix_props.n_collocation)
% instead of randomizing indices of actual meshes, why not introduce a vector that is a random permutation of 1:n_collocation and then here, use
% randvec(global_node_ind) to get the actual collocation pt index that we get vert coords etc for.  This randomizes order of loop over collocation pts but
% keeps meshes in original (presumably somewhat logical) node order.  After matrix solve, can un-randomize order of solution vector so that
% traction, free slip velocity, etc match the original meshes.

index_mapping = struct('local_node2global_node',{struct2cell(index_mapping.local_node2global_node)},...
    'global_node2local_node',{reshape( struct2cell(index_mapping.global_node2local_node) , [] ,1)},...
    'global_node2global_u',index_mapping.global_node2global_u,...
    'global_u2global_node',index_mapping.global_u2global_node);

% Network.link_members = reshape( struct2cell( Network.link_members), [], 1);
Network = struct('n_nodes',Network.n_nodes,'nodes',Network.nodes,'n_links',Network.n_links,'links',Network.links,'E',Network.E,'l_0',Network.l_0,...
    'l',Network.l,'eta',Network.eta,'link_members',{reshape( struct2cell( Network.link_members), [], 1)},'g',Network.g);
% global2local is a n_global_nodes cell array, each cell is a matrix with rows = submesh, local node ind, element ind, position in element
%         index_mapping.global_node2local_node{global_ind}(end+1:end+size(element_inds,1),:) = [repmat([submesh_ind, local_node_ind],size(element_inds,1),1), element_inds, position_inds];
% switch BI_parameters.x0_location
%     case "on_mesh"
%         temp = unique( index_mapping.global_node2local_node{global_node_ind}(:,[1 2]), 'rows' ); % submeshes and local node inds for this global node
%         coll_local = struct('submeshes',temp(:,1),'nodes',temp(:,2));
%     case "off_boundary"
%         coll_local = struct('submeshes',NaN(0,1),'nodes',NaN(0,2));
% end
%     coll_local.submeshes = temp(:,1);
%     coll_local.nodes = temp(:,2);


%     for j = 1:length(Mesh)
%         coll_i_local = find(global_node_ind == Mesh(j).indices.glob.vert,1);
%         if ~isempty(coll_i_local)
%             mesh_ind_coll = j; % index of submesh that current collocation pt belongs to (theoretically arbitrarily chosen between multiple submeshes
%             % that share the coll pt)
%             break
%         end
%         if isempty(coll_i_local) && j == length(Mesh)
%             error('Can''t find global collocation index in any submesh');
%         end
%     end

%     eps2 = BI_parameters.accuracy.eps2(mesh_ind_coll); % eps^2 going with reg. Stokeslets for all collocation points on current submesh
% innertic = tic;

AE = [];   NV = [];     FL = logical([]); %stop parfor complaints about temp variables
r_nodes = []; u_element_nodes = []; RHS_coll0 = []; BI_coll0 = [];
%             local_vert = local_verts(global_node_ind);  replaced by below, if it doesn't work may have to go back to original
%             [~, coll_i_local] = global2local(global_node_ind, Mesh, vert_range, 'vert');


% x_coll = Mesh(coll_local.submeshes(1)).nodes(coll_local.nodes(1),:)';

integrand_constants = struct('eps2',NaN,...
    'x_coll',x0_parameters.x0,...
    'element_nodes',NaN(6,3),...
    'shape_parameters',NaN(3,1),...
    'integral_type',"",...
    'DL_singularity_removal',false,... % doesn't matter for single layer integral, will overwrite later for double layer integral
    'refpoint',NaN(3,1),...
    'field_vel',NaN(6,3)); % either 1 - 6 or 0 for if/where collocation pt is on current element, only matters for DL singularity removal trick with DL integrals
% as with x0_parameters.u, shouldn't matter which instance of this node we look up the coordinates for, should all be the same
% eps is allowed to vary between submeshes so it is defined later in the loop over elements
% 'coll_node_ind', [],...

% r_col = x0_parameters.x0 - BI_parameters.constants.refpoint; % vector from refpoint to collocation pt
% switch BI_parameters.x0_location
%     case "on_mesh"
%         BC_type_coll = mesh_node_parameters.BC_type(global_node_ind);
%         % index for the unknown velocity at the current collocation point,
%         % between 1 - n where n is number of free slip nodes
%         global_u_ind = index_mapping.global_node2global_u( global_node_ind );
%     case "off_mesh"
%         BC_type_coll = NaN;
%         global_u_ind = NaN;
% end



%%temporary matrix for boundry integrals (BI_coll), row is x, y, z components for current collocation point BIE and columns are coeffs for traction components at each vert, to later insert into submatrices
BI_coll = zeros(3,matrix_props.n_cols); % 3 rows for x, y, z components of BIE, columns for x y z components of traction followed by x y z components of unknown u (for free slip coll pts) followed by any unknown kinematic variables (U, Omega, omega)
RHS_coll = zeros(3,matrix_props.n_RHS); % just for current collocation pt, for all flowcases

if BI_parameters.performance.rigid_body_matrix_rotation && firstrun
    BI_coll0 = zeros(3,matrix_props.n_cols);
    RHS_coll0 = zeros(3,matrix_props.n_RHS);
    % else
    %     BI_coll0 = [];  RHS_coll0 = [];
end

% there may be a small penalty to having only 3 rows and many columns (since arrays are stored column-major) but seems unlikely this is significant
% since these arrays aren't *that* big?  Can check where time is spent with Profiler before mexing?
% also, if we flipped the dimensions for these, we'd have to do a transpose later to insert them into the larger arrays, which would be
% expensive, so it's a matter of which penalty is larger...

%A_y_temp = zeros(matrix_props.n_col*3,1); A_z_temp = zeros(matrix_props.n_col*3,1);
%             A_temp = BI_coll; %what actually gets inserted into A - is the same as BI_coll for no-slip BCs but will be different for free-slip



%         RHS_temp0 = zeros(3,matrix_props.n_RHS);
for mesh_ind_elem = 1:length(Mesh)

    
    % if coll pt is a member of the element-containing submesh, then can do tensor rotations since the coll pt and element must be rigidly connected, even if
    % another submesh the coll pt belongs to is deforming
    if BI_parameters.performance.rigid_body_matrix_rotation && ~firstrun && Mesh(mesh_ind_elem).is_mesh_rigid && ismember( mesh_ind_elem , x0_parameters.submeshes )
        % we can simply rotate the tensor contributions coming from the nodes on this element instead of integrating again
        continue
    end
    
        % ignore interaction if the element submesh is not a member of any submeshes the coll pt is on
    if BI_parameters.accuracy.mesh.ignore_interaction && ~ismember( mesh_ind_elem, x0_parameters.submeshes )
        
        %collocation point and element are in different submeshes - skip integrals
        continue
    end
    
    integrand_constants.eps2 = BI_parameters.accuracy.mesh.eps2(mesh_ind_elem); % each submesh has an eps corresponding to all its elements
    % since we're assuming eps is piecewise constant across elements, there can be discontinuities between submeshes, but this should be OK?
    
    for local_elem = 1:Mesh(mesh_ind_elem).n_elements
        %     for elem_i = 1:tot_elems  %elem_i is a global index
        %         elem_i / Mesh.n_elem
        
        
        BI_elem = zeros(3,matrix_props.n_cols); % stores cumulative matrix contributions for integrals over a single element
        RHS_elem = zeros(3,matrix_props.n_RHS); % stores cumulative RHS contributions for integrals over a single element
        
        
        %want mesh_ind_elem below since we are generally interested in current
        %element; node only matters as far as it's distance away when computing reg stokeslet
        %function
        element_node_inds = Mesh(mesh_ind_elem).elements(local_elem,:); % local inds for nodes of current element  1 x 6
        element_nodes = Mesh(mesh_ind_elem).nodes(element_node_inds,:); % coords of nodes of current element   6 x 3
        element_global_node_inds = index_mapping.local_node2global_node{mesh_ind_elem}(element_node_inds); % global inds for nodes of current element  6 x 1
        
        %         element_global_node_inds = Mesh(mesh_ind_elem).indices.glob.vert(element_node_inds); % global inds for verts of current element
        %         submesh_members = {Mesh(mesh_ind_elem).indices.submesh_members{element_node_inds}}; % submeshes each vert of this element belong to
        
        
        BC_types = mesh_node_parameters.BC_type( element_global_node_inds ); % BC for each element node:  1 for no slip, 2 for free slip
        if ~BI_parameters.performance.eliminate_DL(mesh_ind_elem) % include DL integrals over current submesh elements?
            u_element_nodes = ( Mesh(mesh_ind_elem).u(element_node_inds,:,:) ); % known u contributions (body frame velocities for mobility problem, fixed frame velocities
            % resistance problem) at each element node.   6 x 3 x F flowcases
        end
        
        if BI_parameters.problemtype == "mobility"
            r_nodes = element_nodes' - repmat(BI_parameters.constants.refpoint,1,6); % 3 x 6 vectors from refpoint to each vert of current element
        end
        
        unknown_u_global_inds = index_mapping.global_node2global_u(element_global_node_inds); % indices into matrix in case there are unknown fluid u due
        % to free slip BC at some nodes of current element
        
        integrand_constants.element_nodes = element_nodes;
        integrand_constants.shape_parameters = Mesh(mesh_ind_elem).shape_parameters(local_elem,:)';  %[alpha beta gamma]'  3 x 1
        
        switch x0_parameters.location
            case "on_mesh"
                integrand_constants.coll_node_ind = find( global_node_ind == element_global_node_inds);  % needed to appease Coder, which requires this be scalar
%             case "off_mesh"
%                 integrand_constants.coll_node_ind = [];
            otherwise
                integrand_constants.coll_node_ind = [];
        end
        
        %scale abs integration tolerance by element area, since big elements
        %should be allowed more total error than small elements
        abstol_stokeslet =  BI_parameters.accuracy.mesh.integration_tol.stokeslet.abstol * Mesh(mesh_ind_elem).area(local_elem);
        %don't scale reltol, since reltol should already account for area
        %since the integral itself should scale with area
        
        %the important part:  adaptively compute boundary integrals of reg stokeslet function over
        %current element, with respect to current collocation point
        integrand_constants.integral_type = "reg_stokeslet";
        
%         if BI_parameters.performance.debug_mode
            % shatlab hangs at almost no CPU usage with the tic/toc
            % inside adsimp so timing inside the parfor doesn't work even though it
            % compiles without errors
            %[SL_elem, AE, NV, FL, Time] = adsimp( 2, reference_nodes, 54,  BI_parameters.accuracy.integration_tol.traction.maxevals,abstol, BI_parameters.accuracy.integration_tol.traction.reltol ,rule  , rule_constants ,element_nodes,shape_parameters,integrand_constants,'reg_stokeslet', true);
%             [SL_elem, AE, NV, FL] = adsimp( 2, BI_parameters.accuracy.triangle_integration.reference_nodes, 54,  BI_parameters.accuracy.mesh.integration_tol.stokeslet.maxevals,abstol_stokeslet, BI_parameters.accuracy.mesh.integration_tol.stokeslet.reltol ,BI_parameters.accuracy.triangle_integration,integrand_constants);
            %                     SL_elem = ones(54,1);  AE = ones(54,1);  NV = 1;  FL = 1;
%         else %don't bother storing extra outputs
            %                     [SL_elem, ~, ~, ~, ~] = adsimp( 2, reference_nodes, 54,  BI_parameters.accuracy.integration_tol.traction.maxevals,abstol, BI_parameters.accuracy.integration_tol.traction.reltol ,rule  , rule_constants ,element_nodes,shape_parameters,integrand_constants,'reg_stokeslet', false);
            
            [SL_elem, ~, ~, ~] = adsimp( 2, BI_parameters.accuracy.triangle_integration.reference_nodes, 54,  BI_parameters.accuracy.mesh.integration_tol.stokeslet.maxevals,abstol_stokeslet, BI_parameters.accuracy.mesh.integration_tol.stokeslet.reltol ,BI_parameters.accuracy.triangle_integration,integrand_constants);
            %                 SL_elem = ones(54,1);  AE = ones(54,1);  NV = 1;  FL = 1;
            SL_elem = SL_elem * 1/2 / BI_parameters.constants.mu;  % assumes that solid angle (2*pi for a smooth region) multiplies u_c on left side of BIE
            
%         end
        %               out =  hS*[S(1,1)*phi; S(1,2)*phi; S(1,3)*phi; ...
        %             S(2,1)*phi; S(2,2)*phi; S(2,3)*phi; ...
        %             S(3,1)*phi; S(3,2)*phi; S(3,3)*phi];
        % integrals of phiS
        submat = reshape(SL_elem,18,3)'; % 6 nodes, 3 x y z components of traction per node = 3 rows for x y z components of BIE,
        % 18 columns for x, y, z traction components and 6 phi
        
        %entries in matrix are total boundary integrals (BI) so can do a cumulative sum
        inds = ( repmat( 3 * (element_global_node_inds - 1),1,3) + repmat([1 2 3],6,1) ); % indices of matrix columns corresponding to x,y,z traction components at global nodes of
        % current element   6 x 3 each column are the matrix column inds of the 6 nodes, and we have 3 columns for x, y, z components of traction
        
        
        BI_elem(:,inds(:)) = BI_elem(:,inds(:)) + submat;
        %                 BI_elem(2,col_inds) = BI_elem(2,col_inds) + submat(:,2)';
        %                 BI_elem(3,col_inds) = BI_elem(3,col_inds) + submat(:,3)';
        
        %             if BI_parameters.performance.rigid_body_matrix_rotation && firstrun && Mesh(mesh_ind_elem).rigid && mesh_ind_elem == mesh_ind_coll
        %                 BI0(:,col_inds) = BI0(:,col_inds) + submat';  % this version of BI_coll will only have constributions from integrals for elements and coll
        %                 % pts on the same rigid body
        %             end
        
        
        if ~BI_parameters.performance.eliminate_DL(mesh_ind_elem)  % compute DL integrals over this element
            integrand_constants.integral_type = "reg_stresslet";
            integrand_constants.DL_singularity_removal = BI_parameters.performance.DL_singularity_removal(mesh_ind_elem);
            abstol_stresslet =  BI_parameters.accuracy.mesh.integration_tol.stresslet.phiTn.abstol * Mesh(mesh_ind_elem).area(local_elem);
            % integrating phiTn regardless of no slip or free slip, but need to handle all flowcases for no slip
            [DL_elem, ~, ~, ~] = adsimp( 2, BI_parameters.accuracy.triangle_integration.reference_nodes, 54,  BI_parameters.accuracy.mesh.integration_tol.stresslet.phiTn.maxevals, abstol_stresslet, BI_parameters.accuracy.mesh.integration_tol.stresslet.phiTn.reltol ,BI_parameters.accuracy.triangle_integration, integrand_constants);
            DL_elem = DL_elem * 1/2;  %  assumes that solid angle (2*pi for a smooth region) multiplies u_c on left side of BIE
            
            switch BI_parameters.problemtype
                case "resistance"
                    
                    %                                 BI_parameters.BC_type(mesh_ind_elem) % only matters for inserting values into A_BIE or RHS, either way we're integrating phiTn
                    % mobility:  for elements that are shared between no slip, free slip vertices, can simply compute all variations of the integrand (e.g.
                    % Tn, rTn, phiTn) and just only use the ones needed for each vertex of the element
                    
                    
                    
                    % each Mesh(i).u is Nx3xM where M is length(flowcases), N is #verts
                    %                             u_colls = zeros(matrix_props.n_col, matrix_props.n_unknown_u*3 ,3); % extra columns for coeffs corresponding to unknown u components in DL integrals
                    %                             u_rows = zeros(matrix_props.n_unknown_u*3 , matrix_props.n_cols);
                    
                    
                    
                    for vi = 1:6
                        % phiTn for phi of current vert
                        phiTn = [ DL_elem(vi)  DL_elem(6+vi) DL_elem(12+vi) ;...
                            DL_elem(18+vi)  DL_elem(24+vi)  DL_elem(30+vi) ;...
                            DL_elem(36+vi) DL_elem(42+vi)  DL_elem(48+vi) ;  ];
                        
                        switch BC_types(vi)
                            
                            case 1 % no slip
                                
                                % flowcases, e.g. forced x y z translation and x y z rotation to calculate diffusivities
                                if BI_parameters.performance.DL_singularity_removal(mesh_ind_elem) && x0_parameters.BC_type == 1
                                    % no slip, x0_parameters.u is prescribed and known
                                    vel_mod = reshape( squeeze(u_element_nodes(vi,:,:)),3,[]) - x0_parameters.u;  %u_element_nodes is 6 x 3 x flowcases, x0_parameters.u is 3 x flowcases
                                    % free slip, x0_parameters.u is unknown
                                else % either not doing trick, or doing trick but x0_parameters.u is unknown
                                    vel_mod = reshape( squeeze(u_element_nodes(vi,:,:)),3,[]);
                                end
                                % velmod is 3 x flowcases
                                % in any case, add integral corresponding to known u_v to RHS
                                RHS_elem = RHS_elem - phiTn * vel_mod; % RHS_elem is 3 x flowcases
                                % minus since we're moving this term from LHS to RHS
                                %                                     if BI_parameters.performance.rigid_body_matrix_rotation && firstrun && Mesh(mesh_ind_elem).rigid && mesh_ind_elem == mesh_ind_coll
                                %                                         RHS_temp0 = RHS_temp0 - phiTn * vel_mod;   % this version of BI_coll will only have constributions from integrals for elements and coll
                                %                                         % pts on the same rigid body
                                %                                     end
                                
                            case 2 % free slip
                                % in case of no trick, these coeffs go with u for this vert (across columns for each u component, down rows for each BIE
                                % component)
                                inds = matrix_props.n_collocation*3 + 3 * (unknown_u_global_inds(vi) - 1) + [1 2 3]; % corresponding to inds for 3 components of unknown u
                                % matrix columns are organized with unknown u following unknown traction
                                BI_elem(:,inds ) = BI_elem(:,inds ) + phiTn;
                                
                                % in case of trick, these coeffs still go with u for this vert, AND -coeffs go with x0_parameters.u:
                                if BI_parameters.performance.DL_singularity_removal(mesh_ind_elem) && x0_parameters.BC_type == 1
                                    % if x0_parameters.u is no-slip, then add sum(coeffs .* x0_parameters.u , 2) to RHS.
                                    % flowcases, e.g. forced x y z translation and x y z rotation to calculate diffusivities
                                    RHS_elem = RHS_elem + phiTn * x0_parameters.u; % minus sign from moving from LHS to RHS negated since we want to add integral of -x0_parameters.u
                                    
                                end
                                
                        end % no slip vs free slip
                        
                        if BI_parameters.performance.DL_singularity_removal(mesh_ind_elem) && x0_parameters.BC_type == 2
                            % if x0_parameters.u is free-slip, also add integral corresponding to unknown u_c to coeffs in matrix
                            inds = matrix_props.n_collocation*3 + 3 * (x0_parameters.global_u_ind - 1) + [1 2 3]; % corresponding to inds for 3 components of unknown u
                            % matrix columns are organized with unknown u following unknown traction
                            % minus sign since we were integrating (u_v - u_c)*Tn
                            BI_elem(:,inds ) = BI_elem(:,inds) - phiTn;
                        end
                        
                    end % element_nodes of current element
                    
                    
                    
                case "mobility"
                    
                    
                    %mobility, noslip   same matrix size    add DL terms into coeffs for U, Omega, and export DL terms for u_r to RHS
                    % integrate Tn for U coeffs (9 values for 3 U components, 3 BIE components), rTn for Omega coeffs (9 values), phiTn to overall integrate uTn
                    % where u = known, interped u_r on flagella (6 x 9 = 54 values, combined with known u_r on element verts to yield 3 values for 3 components of BIE, exported to RHS)
                    % in trick, Tn is set to zeros since int(U - U_c)Tn will be zero everywhere
                    % in trick, we overall int( Omega x (r - r_c) Tn ) so that integrand goes to zero at collocation vertex
                    % in trick, we manually zero out phiTn when an element vertex is the same as the current collocation pt
                    %                                 maxevals_temp = [ repmat(BI_parameters.accuracy.integration_tol.stresslet.Tn.maxevals,9,1); repmat(BI_parameters.accuracy.integration_tol.stresslet.rTn.maxevals,9,1); repmat(BI_parameters.accuracy.integration_tol.stresslet.phiTn.maxevals,54,1); ];
                    %                                 abstol_temp = [ repmat(abstol.stresslet.Tn,9,1); repmat(abstol.stresslet.rTn,9,1);  repmat(abstol.stresslet.phiTn,54,1); ];
                    %                                 reltol_temp = [ repmat(BI_parameters.accuracy.integration_tol.stresslet.Tn.reltol,9,1); repmat(BI_parameters.accuracy.integration_tol.stresslet.rTn.reltol,9,1); repmat(BI_parameters.accuracy.integration_tol.stresslet.phiTn.reltol,54,1); ];
                    
                    %                                 [DL_elem, ~, ~, ~] = adsimp( 2, reference_nodes, 54,  maxevals_temp, abstol_temp,reltol_temp ,rule  , rule_constants ,element_nodes,shape_parameters,integrand_constants,'reg_stresslet',{'Tn','rTn','phiTn'},BI_parameters.performance.DL_singularity_removal);
                    %                                 Tn = reshape(DL_elem(1:9),3,3);  rTn = reshape(DL_elem(10:18),3,3);  phiTn = DL_elem(19:72);
                    
                    for vi = 1:6
                        % phiTn for phi of current vert
                        phiTn = [ DL_elem(vi)  DL_elem(6+vi) DL_elem(12+vi) ;...
                            DL_elem(18+vi)  DL_elem(24+vi)  DL_elem(30+vi) ;...
                            DL_elem(36+vi) DL_elem(42+vi)  DL_elem(48+vi) ;  ];
                        
                        switch BC_types(vi)
                            
                            case 1 % no slip
                                
                                % keeping the ability to handle multiple "flowcases" here for giggles - i.e. if one wanted to try
                                % relative surface velocities, perhaps for a squirmer? (geometry would need to be the same over time
                                % for each case, which is kind of limiting....)
                                if BI_parameters.performance.DL_singularity_removal(mesh_ind_elem) && x0_parameters.BC_type == 1
                                    % no slip, x0_parameters.u is prescribed and known
                                    vel_mod = reshape( squeeze(u_element_nodes(vi,:,:)),3,[]) - x0_parameters.u;  %u_element_nodes is 6 x 3 x flowcases, x0_parameters.u is 3 x flowcases
                                    r_mod = r_nodes(:,vi) - x0_parameters.r;
                                    U_factor = 0; % int(U_v - U_c)Tn is always zero, even when current node is not same as current coll pt, since U is
                                    % a constant everywhere
                                    % free slip, x0_parameters.u is unknown
                                else % either not doing trick, or doing trick but x0_parameters.u is unknown
                                    vel_mod = reshape( squeeze(u_element_nodes(vi,:,:)),3,[]);
                                    r_mod = r_nodes(:,vi);
                                    U_factor = 1;
                                end
                                % vel_mod is 3 x flowcases
                                
                                % rphiTn are the coeffs going with (Omega x r)phiTn   3 x 3
                                rphiTn = [r_mod(2)*phiTn(:,3) - r_mod(3)*phiTn(:,2) , r_mod(3)*phiTn(:,1) - r_mod(1)*phiTn(:,3) , r_mod(1)*phiTn(:,2) - r_mod(2)*phiTn(:,1) ];
                                % may have to repmat each r into 3 rows for Coder
                                
                                % these entries hold whether or not we do the trick since the trick is taken care of when we compute Tn, rTn and with
                                % U_factor and r_mod
                                inds = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + (1:6);
                                BI_elem(:,inds ) = BI_elem(:,inds ) + [repmat(U_factor,3,3).*phiTn  rphiTn]; % coeffs for U and Omega
                                % need U_factor here to zero out the U phiTn term with trick even when current node is not same as current coll pt, since U phiTn
                                % term always cancels out when doing trick
                                
                                if BI_parameters.rotating_flagellum && strcmp(BI_parameters.Tail.motorBC,'torque') && mesh_ind_elem == BI_parameters.Tail.submesh_index
                                    % (only accounting for one tail for now, need to generalize for multiflagellated)
                                    % if element node is on bacterial tail and no slip, with unknown omega
                                    ind = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 7; % coeff for unknown omega for tail rotation relative to body
                                    BI_elem(:,ind ) = BI_elem(:,ind ) - rphiTn * BI_parameters.Tail.motor_orientation; % see "double layer" work
                                    % negative sign due to how we define the contribution of relative tail rotation (see Schuech et al PNAS SI)
                                end
                                
                                
                                % in any case, add integral corresponding to known u_v to RHS
                                RHS_elem = RHS_elem - phiTn * vel_mod; % RHS_elem is 3 x flowcases
                                % minus since we're moving this term from LHS to RHS
                                
                            case 2 % free slip
                                % in case of no trick, these coeffs go with u for this vert (across columns for each u component, down rows for each BIE
                                % component)
                                inds = matrix_props.n_collocation*3 + 3 * (unknown_u_global_inds(vi) - 1) + [1 2 3]; % corresponding to inds for 3 components of unknown u
                                % matrix columns are organized with unknown u following unknown traction
                                BI_elem(:,inds ) = BI_elem(:,inds ) + phiTn;
                                
                                % in case of trick, these coeffs still go with u for this vert as included directly above, AND -coeffs also go with x0_parameters.u:
                                if BI_parameters.performance.DL_singularity_removal(mesh_ind_elem) && x0_parameters.BC_type == 1
                                    % need to include -int(U phiTn) - int(Omega x r phiTn) for collocation pt
                                    % rphiTn are the coeffs going with (Omega x r)phiTn   3 x 3
                                    rphiTn = [x0_parameters.r(2)*phiTn(:,3) - x0_parameters.r(3)*phiTn(:,2) , x0_parameters.r(3)*phiTn(:,1) - x0_parameters.r(1)*phiTn(:,3) , x0_parameters.r(1)*phiTn(:,2) - x0_parameters.r(2)*phiTn(:,1) ];
                                    % may have to repmat each r into 3 rows for Coder
                                    
                                    inds = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + (1:6);
                                    BI_elem(:,inds ) = BI_elem(:,inds ) - [phiTn  rphiTn]; % coeffs for U and Omega
                                    
                                    if BI_parameters.rotating_flagellum && BI_parameters.Tail.motorBC == "torque" && ismember( BI_parameters.Tail.submesh_index , x0_parameters.submeshes )
                                        ind = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 7; % coeff for unknown omega for tail rotation relative to body
                                        BI_elem(:,ind ) = BI_elem(:,ind ) + rphiTn * BI_parameters.Tail.motor_orientation; % see "double layer" work
                                        % plus:  normally would be minus due to how tail rotation is accounted for, see Schuech PNAS 2019, but we are
                                        % subtracting x0_parameters.u via the trick so it becomes plus
                                    end
                                    
                                    
                                    % also add int(u_r phiTn) at collocation pt to RHS.
                                    % flowcases, e.g. forced x y z translation and x y z rotation to calculate diffusivities
                                    RHS_elem = RHS_elem + phiTn * x0_parameters.u; % minus sign from moving from LHS to RHS negated since we want to add integral of -x0_parameters.u
                                    
                                end
                                
                        end % no slip vs free slip
                        
                        if BI_parameters.performance.DL_singularity_removal(mesh_ind_elem) && x0_parameters.BC_type == 2
                            % if x0_parameters.u is free-slip, also add integral corresponding to unknown u_c to coeffs in matrix
                            inds = matrix_props.n_collocation*3 + 3 * (x0_parameters.global_u_ind - 1) + [1 2 3]; % corresponding to inds for 3 components of unknown u
                            % matrix columns are organized with unknown u following unknown traction
                            % minus sign since we were integrating (u_v - u_c)*Tn
                            BI_elem(:,inds ) = BI_elem(:,inds) - phiTn;
                        end
                        
                        
                    end % element_nodes of current element
                    
                    
                    
            end % resistance problem vs mobility problem
            
        end % include DL
        
        
        BI_coll = BI_coll + BI_elem;
        RHS_coll = RHS_coll + RHS_elem;
        
        if BI_parameters.performance.rigid_body_matrix_rotation && firstrun && Mesh(mesh_ind_elem).is_mesh_rigid && ismember(mesh_ind_elem , x0_parameters.submeshes)
            
            BI_coll0 = BI_coll0 + BI_elem;
            RHS_coll0 = RHS_coll0 + RHS_elem;   % this version of BI_coll will only have contributions from integrals for elements and coll
            % pts on the same rigid body
        end
        
        
        
    end  % elements
    
end % submeshes

% now account for "boundary integrals" over 0D network nodes (just a summation of force*S at each node, no integrals of traction needed)

temp = zeros(3,1);
for network_node_ind = 1:Network.n_nodes
    
    temp = temp + calcS(Network.nodes(network_node_ind,:)',x0_parameters.x0,BI_parameters.accuracy.network.eps2) * ( Network.g(network_node_ind,:) + Repulsion.F(network_node_ind,:) )' ;

    %% make sure this is zero for no links
end

RHS_coll = RHS_coll - 1/BI_parameters.constants.mu * 1/2 * temp; % moving sum over Stokeslets from LHS to RHS, so subtract

% will never be able to do any recycling via RHS_coll0 for network nodes since they will always move relative to all other nodes?