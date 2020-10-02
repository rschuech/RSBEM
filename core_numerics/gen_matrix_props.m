function [matrix_props] = gen_matrix_props(input,Mesh,node_parameters)

matrix_props.n_collocation = length(node_parameters.node_type); 

matrix_props.n_unknown_u = sum(node_parameters.BC_type == 2);

switch input.problemtype
    case 'mobility'
        if input.rotating_flagellum
                switch input.Tail.motorBC
                    case 'freq'
                        n_unknown_kinematics = 6; % U, Omega with 3 components each
                    case 'torque'
                        n_unknown_kinematics = 7; % U, Omega with 3 components each, plus scalar omega for tail rotation rate
                end
        else
                n_unknown_kinematics = 6; % U, Omega with 3 components each
        end
    case 'resistance'
        n_unknown_kinematics = 0; % all velocities are prescribed
        
end


matrix_props.n_rows = matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + n_unknown_kinematics;
% 3 components of traction for each collocation pt, 3 components of each unknown u due to free slip BC, an equation for each unknown kinematic variable


matrix_props.n_cols = matrix_props.n_rows; % A should always be square for the forseeable future
% matrix_props.Col_inds = save_Col_inds_mex(Mesh,input.performance.nthreads);  % not using mexed for now, error due to # fields in Mesh now different

matrix_props.n_RHS = size(Mesh(1).u,3); % number of RHS to solve system with

% %store matrix columns (i.e., x y z traction) going with each node
%     tot_nodes = Mesh(end).indices.glob.unq_bounds.node(2 ); % total # global verts
%     temp = repmat(1:tot_nodes,3,1); % [1:N; 1:N; 1:N] for x y z components for all verts
%     vertind_cols = temp(:)';  %defined in relation to global indices  [x y z x y z x y z ...] or [1 1 1 2 2 2 3 3 3 ... N N N] where N is final global vert
