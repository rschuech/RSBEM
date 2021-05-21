function [A, A_force, A_torque, RHS, A_motor_torque] = matrix_assembly_mex_wrapper(Mesh,Network, Repulsion, matrix_props,index_mapping,mesh_node_parameters,assembly_input, t)
%wraps mexed function and assembles full A in native Matlab to avoid
%ridiculous limit for matrix size in Coder
%%
%make sure we can even initialize the full damn thing before going further
A = zeros(matrix_props.n_rows, matrix_props.n_cols); % entire matrix 
RHS = zeros(matrix_props.n_rows, matrix_props.n_RHS);

% switch assembly_input.problemtype
%     case "resistance"
% index_mapping2 = index_mapping;
% index_mapping2.local_node2global_node = cell2struct(index_mapping.local_node2global_node,'indices',1);
% index_mapping2.global_node2local_node = cell2struct(index_mapping.global_node2local_node,'indices',length(index_mapping.global_node2local_node));

% rng(0);
% rand_row_inds = randperm(matrix_props.n_collocation)';
 rand_row_inds = 1:matrix_props.n_collocation;
rand_row_inds = [rand_row_inds; rand_row_inds + matrix_props.n_collocation; rand_row_inds + 2*matrix_props.n_collocation; ];



% A and RHS are now initialized inside matrix_assembly_mex but this can be easily done out here, then send them in and get them back out, modified 
% in-place?
        [A(rand_row_inds,:), RHS(rand_row_inds,:)] = matrix_assembly_mex(Mesh,Network, Repulsion, matrix_props,index_mapping,mesh_node_parameters,assembly_input);
        % pretty sure this unfortunately has A_BIE and A separately in memory but not sure how to avoid the two copies
        % because of parfor rules, we have to assemble the big BIE matrix block as one variable that can be reshaped between 3D and 2D - the reshaping does
        % not work if we include the extra rows and columns of full A
        % so somewhere in the code, we have to have the BIE matrix block as one variable and the full matrix as another variable, no easy way around it?
        % or, maybe can do a permute-type operation instead of reshape on the BIE block but that would create a temporary copy anyway....
        
%     case "mobility"
%         [A_BIE_1,  A_BIE_2,  A_BIE_3,  A_BIE_4,  RHS(1:matrix_props.n_collocation*3,:), debug_info] = matrix_assembly_mex(Mesh,matrix_props,assembly_input);
% end

% if ~isempty(A_BIE_2), inds_2 = [1:3]; else, inds_2 = [1 1 1]; end % one is allowed to ask for the first page of an empty matrix, which yields an empty matrix
% if ~isempty(A_BIE_3), inds_3 = [1:3]; else, inds_3 = [1 1 1]; end
% if ~isempty(A_BIE_4), inds_4 = [1:3]; else, inds_4 = [1 1 1]; end

% A(1:matrix_props.n_collocation*3, :) = [A_BIE_1 A_BIE_2 A_BIE_3 A_BIE_4]; % combine parts of A here outside mex function as workaround to dumb Coder limitation on array size
% clear A_BIE* % does this help or hurt runtime?
%   A_BIE_1 = zeros(matrix_props.n_collocation,3,matrix_props.n_cols);

% A(1:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = [A_1(:,:,1) A_2(:,:,inds_2(1)) A_3(:,:,inds_3(1)) A_4(:,:,inds_4(1))];  %A is split up over columns
% A(2:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = [A_1(:,:,2) A_2(:,:,inds_2(2)) A_3(:,:,inds_3(2)) A_4(:,:,inds_4(2))];
% A(3:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = [A_1(:,:,3) A_2(:,:,inds_2(3)) A_3(:,:,inds_3(3)) A_4(:,:,inds_4(3))];

%   A = assembly_input.constants.multfactor*A;  %multiply by 1/8/pi as per eq. 3.4, DJ Smith "A boundary element regularized..." paper
  
  
%   if matrix_props.n_unknown_u > 0 % we have free-slip nodes
      rows = matrix_props.n_collocation*3+1  :  (matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3); % row indices for additional free-slip constraints
      [A(rows, :), RHS(rows, :)] = free_slip_constraints(Mesh, matrix_props,index_mapping,mesh_node_parameters,assembly_input);
%   end
  
%compute A_force and A_torque, either for direct output and
%later use in resistance problem case, or immediate insertion into A matrix for
%free-swimming case
%tic

% [~, A_force]  = compute_force_integral(Mesh, matrix_props, assembly_input); %parallelized   A_force is separated via 3rd dimension for each submesh component
% [~, A_torque] = compute_torque_integral(Mesh, matrix_props, assembly_input); %parallelized  A_torque is separated via 3rd dimension for each submesh component
%for all cases but particularly the resistance problem case, A_force and A_torque can
%be multiplied by f to give the total force and torque on the object

% A_force, A_torque are separated via 3rd dimension for each submesh component.  Function is parallelized via parfor over elements of each submesh.
[A_force, A_torque] = compute_force_torque_matrices(Mesh, matrix_props, index_mapping,assembly_input);


if assembly_input.rotating_flagellum && ismember("Tail",[Mesh.name]) % would have to generalize for multiple rotating flagella or tails with multiple submeshes
    tail_torque = A_torque(:,:,[Mesh.name] == "Tail"); % torque between body and tail - should equal the motor torque
    motor_orientation = Mesh([Mesh.name] == "Tail").orientation(:,1);  % updated to allow for tails attached at an angle - motor rotation axis should always be the centerline axis of the tail
    % tail_torque = 3 x N, motor_orientation = 1 x 3, A_motor_torque = 1 x N
%     A_motor_torque = motor_orientation(1)*tail_torque(1,:) + motor_orientation(2)*tail_torque(2,:) + motor_orientation(3)*tail_torque(3,:);
    A_motor_torque = motor_orientation' * tail_torque;
    % multiplied by traction, this gives motor torque:   A_motor_torque * f = motor torque
else
A_motor_torque = [];
end


% A_constraints and RHS_contraints will contain all equations required besides BIEs, e.g. free-slip BC constraints at relevant nodes as well as
% free-swimming force- and torque-free constraints, and for problems with a rotating flagellum with specified torque, the constraint setting that motor
% torque

%%  additional free swimming equations

if assembly_input.problemtype == "mobility"
    rows = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 1 : matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 7; 
    % assume 3 U components, 3 Omega components, and omega even if we aren't solving for that
    A(rows(1:3),  1:matrix_props.n_collocation*3) = sum(A_force,3); %sum over submeshes for integrals over entire mesh
    A(rows(4:6),  1:matrix_props.n_collocation*3) = sum(A_torque,3);
    % no need to explicitly do anything for RHS here for a mobility problem since it's initialized to zeros already
    % except now we do, with external repulsion forces
    
    filter = ~isnan(Repulsion.mesh_index);
    
    RHS(rows(1:3),:) = repmat( - sum(Repulsion.F(filter,:), 1)' , 1, size(RHS,2));
    [refpoint] = get_rotational_references(Mesh, assembly_input);
    
    RHS(rows(4:6),:) = repmat( - sum( cross( (Repulsion.x(filter,:)  - refpoint') , Repulsion.F(filter,:) ), 1 )', 1, size(RHS,2));
    
    % sum of torques on tail  = specified motor torque (only need for
    % last equation of torque specified BC
    
    %do this no matter what and output row of A below for
    %rotationrate BC, but don't actually add row to A
    
    if assembly_input.rotating_flagellum && ismember("Tail",[Mesh.name]) && assembly_input.Tail.motorBC == "torque"
        tail_ind = find("Tail" == [Mesh.name]);
        A(rows(7),  1:matrix_props.n_collocation*3) = A_motor_torque;
        
        filter = Repulsion.mesh_index == tail_ind;
        %        motor_orientation * (- sum( cross( (Repulsion.x(filter,:)  - refpoint') , Repulsion.F(filter,:) ) )'); % repulsive torque on tail (vector)
        RHS(rows(7),:) = assembly_input.Tail.motor_torque(t) + motor_orientation' * (- sum( cross( (Repulsion.x(filter,:)  - refpoint') , Repulsion.F(filter,:) ), 1 )');
%      RHS(rows(7),:) = assembly_input.Tail.motor_torque(t);
    
    end
    
    
    
end  %mobility problem





% switch assembly_input.BC_type
%     case 1
%         A = assembly_input.constants.multfactor*A;  %multiply by 1/8/pi as per eq. 3.4, DJ Smith "A boundary element regularized..." paper
%     case 2
%         A(1:3:end,:) = A(1:3:end,:) * assembly_input.constants.multfactor;
% end


% switch assembly_input.problemtype
%     case 'freeswim'
%         switch assembly_input.bugtype
%             case 'bacteria'
%                 switch assembly_input.tail.motorBC
%                     case 'freq'
%                         A(matrix_props.n_col*3 + 1:matrix_props.n_col*3 + 6, :) = A_freeswim_rows;
%                         A(:, matrix_props.n_col*3 + 1 : matrix_props.n_col*3 + 6) = A_freeswim_cols;
%                     case 'torque'
%                         A(matrix_props.n_col*3 + 1:matrix_props.n_col*3 + 7, :) = A_freeswim_rows;
%                         A(:, matrix_props.n_col*3 + 1 : matrix_props.n_col*3 + 7) = A_freeswim_cols;
%                 end
%             case {'dino', 'sheet'}
%                 A(matrix_props.n_col*3 + 1:matrix_props.n_col*3 + 6, :) = A_freeswim_rows;
%                 A(:, matrix_props.n_col*3 + 1 : matrix_props.n_col*3 + 6) = A_freeswim_cols;
%         end
%         
%         
% end
