function [u_field] = field_velocity(Mesh,Network, Repulsion, x_field, solution, matrix_props,index_mapping,mesh_node_parameters,assembly_input)

[refpoint, motor_orientation] = get_rotational_references(Mesh, assembly_input);


% index_mapping = struct('local_node2global_node',{struct2cell(index_mapping.local_node2global_node)},...
%     'global_node2local_node',{reshape( struct2cell(index_mapping.global_node2local_node) , [] ,1)},...
%     'global_node2global_u',index_mapping.global_node2global_u,...
%     'global_u2global_node',index_mapping.global_u2global_node);



coder.extrinsic('tic');
coder.extrinsic('toc');
coder.extrinsic('num2str');

index_mapping2 = index_mapping;
index_mapping2.local_node2global_node = cell2struct(index_mapping.local_node2global_node,'indices',1);
index_mapping2.global_node2local_node = cell2struct(index_mapping.global_node2local_node,'indices',length(index_mapping.global_node2local_node));

Network2 = Network;
Network2.link_members = cell2struct(Network2.link_members,'indices',length(Network2.link_members));


tic
%#codegen

%numels_max = double(intmax('int32'));  %dumb Coder limit on numel of any array

  
% BI_parameters = assembly_input;
% BI_parameters.Tail = rmfield(BI_parameters.Tail,'motor_torque'); % don't need this for BIEs, only comes in on RHS of additional constraint eqs later
% % and we want to generalize it as a function of time
% BI_parameters = rmfield(BI_parameters,'repulsion');

BI_parameters.performance.eliminate_DL = assembly_input.performance.eliminate_DL;
BI_parameters.performance.DL_singularity_removal = assembly_input.performance.DL_singularity_removal;
BI_parameters.performance.rigid_body_matrix_rotation = assembly_input.performance.rigid_body_matrix_rotation;
BI_parameters.performance.verbose = assembly_input.performance.verbose;
BI_parameters.problemtype = assembly_input.problemtype;
BI_parameters.rotating_flagellum = assembly_input.rotating_flagellum;
BI_parameters.Tail.motorBC = assembly_input.Tail.motorBC;
BI_parameters.Tail.submesh_index = assembly_input.Tail.submesh_index;
BI_parameters.Tail.motor_orientation = assembly_input.Tail.motor_orientation;
% BI_parameters.accuracy = assembly_input.accuracy;
BI_parameters.accuracy.mesh = assembly_input.accuracy.mesh;
BI_parameters.accuracy.network = assembly_input.accuracy.network;
BI_parameters.accuracy.triangle_integration = assembly_input.accuracy.triangle_integration;
BI_parameters.constants = assembly_input.constants;



x0_location = "off_mesh";
BI_parameters.constants.refpoint = refpoint;
BI_parameters.Tail.motor_orientation = motor_orientation;

% coder.varsize('reference_nodes',[2 3 Inf],[true true true]); %for mex rules
u_field = NaN(size(x_field,1),3);

parfor (x_field_ind = 1:size(x_field,1) , assembly_input.performance.nthreads)  %rows of A_BIE in sets of 3, i.e. x y z components of BIE for each node
    % this loops over collocation points
%         for x_field_ind = 1:size(x_field,1)   %randperm(matrix_props.n_collocation)
    % instead of randomizing indices of actual meshes, why not introduce a vector that is a random permutation of 1:n_collocation and then here, use
    % randvec(x_ind) to get the actual collocation pt index that we get vert coords etc for.  This randomizes order of loop over collocation pts but
    % keeps meshes in original (presumably somewhat logical) node order.  After matrix solve, can un-randomize order of solution vector so that
    % traction, free slip velocity, etc match the original meshes.
    
    % need to use below again in separate code that calcs u at network nodes
%     switch BI_parameters.x0_location
%         case "on_mesh"
%             temp = unique( index_mapping.global_node2local_node{x_ind}(:,[1 2]), 'rows' ); % submeshes and local node inds for this global node
%             coll_local = struct('submeshes',temp(:,1),'nodes',temp(:,2));
%         case "off_boundary"
%             coll_local = struct('submeshes',NaN(0,1),'nodes',NaN(0,2));
%     end
    
%       temp = unique( index_mapping.global_node2local_node{x_ind}(:,[1 2]), 'rows' ); % submeshes and local node inds for this global node
            x0_parameters = struct('x0',x_field(x_field_ind,:)',  'submeshes',NaN,      'nodes',NaN);
    x0_parameters.location = x0_location;
            
     x0_parameters.r = x_field(x_field_ind,:)' - BI_parameters.constants.refpoint; % vector from refpoint to collocation pt
            
            
%             x0_parameters.r = Mesh(temp(1,1)).nodes(temp(1,2),:)' - BI_parameters.constants.refpoint; % vector from refpoint to collocation pt
% switch x0_location
%     case "on_mesh"
%         x0_parameters.BC_type = mesh_node_parameters.BC_type(x_ind);
%         % index for the unknown velocity at the current collocation point,
%         % between 1 - n where n is number of free slip nodes
%         x0_parameters.global_u_ind = index_mapping.global_node2global_u( x_ind );
%         x0_parameters.u = reshape(  squeeze(Mesh(temp(1,1)).u(temp(1,2),:,:)), 3 , []  ); % needed if including DL  3 x flowcases
%
%     case "off_mesh"
x0_parameters.BC_type = NaN; % actually no slip but this should stop code from ever trying singularity removal trick regardless of singularity removal trick flag
x0_parameters.global_u_ind = NaN;
x0_parameters.u = NaN(3,1); % should only ever be used in conjunction with singularity removal trick, which we can't use if trying to directly compute u at x0 (?)
% end


% reshape needed in case there's only one flowcase, to make sure
% arbitrarily look up u_coll for first submesh/local node this global node appears in - u_coll should be the same for all appearances


    firstrun = true;
   [BI_coll, RHS_coll] = boundary_integrals_mexed(x0_parameters, NaN, Mesh, mesh_node_parameters, Network2, Repulsion, matrix_props, index_mapping2, BI_parameters, firstrun );
% RHS_coll will be zero if eliminating the DL, since that's the only contribution to the RHS besides from the original u-side of the BIE, which is not incorporated
% here since we're solving for it when we calculate field velocities
  

% diffs = A([x_field_ind, matrix_props.n_collocation + x_field_ind, matrix_props.n_collocation*2 + x_field_ind],1:matrix_props.n_collocation*3) - BI_coll(:,1:matrix_props.n_collocation*3);
% diffs = A([x_field_ind, matrix_props.n_collocation + x_field_ind, matrix_props.n_collocation*2 + x_field_ind],:) - BI_coll(:,:);


% if max(abs(diffs(:))) > 1E-12
%     stopa
% end

u_field(x_field_ind,:) = (BI_coll*solution - RHS_coll) / (4*pi);  % hard coding alpha = 4*pi at the field point here, i.e. assuming it is not on a solid surface
    
% shat = BI_coll(:,end-6:end);
% 
%  if any(shat(:) ~= 0)
%      stoapa
%  end


% if any(RHS_coll) ~= 0
%     stopa
% end


end   %parfor over collocation points %%%%%%%%%%%%%%%%%%%%%%%%%    parfor      %%%%%%%%%%%%%%%%%%%%%%%%%

 
