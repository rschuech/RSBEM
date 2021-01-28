function [dydt, Mesh, Network, Repulsion] = derivatives(t,y, Mesh, Network, mesh_node_parameters, index_mapping, matrix_props, assembly_input)

%% Update Mesh and Network based on y
[Mesh, Network, mesh_node_parameters, y_swimmer, r] = update_Mesh_Network(Mesh, Network, t, y, mesh_node_parameters,index_mapping, assembly_input);
 
% is_broken = Network.E == 0; % all broken links, both that were already broken and that just broke
% [Mesh_bodyframe, mesh_node_parameters_bodyframe, A_1_bodyframe] = rotateMesh(Mesh,Mesh(2).refpoints(:,1),mesh_node_parameters,index_mapping, ...
%     2, y(end),Mesh(2).orientation);
%%
Repulsion = calc_repulsive_forces(Mesh, Network, index_mapping, assembly_input.repulsion);
% [min_dist, Repulsion.x, xi_eta, Repulsion.mesh_index, element_index, is_inside] = dist2mesh(Network.nodes, Mesh, assembly_input.repulsion.mindist2); 
% too_close = min_dist < repulsion.d;
% vec = Network.nodes(too_close,:) - Repulsion.x(too_close,:);
% dir = vec ./ sqrt(  sum( ( vec ).^2 , 2 ) );
% Repulsion.F = zeros(Network.n_nodes,3);
% % repulsion forces on network nodes
% Repulsion.F(too_close,:) = repulsion.g * ( exp(-min_dist(too_close)/repulsion.d) - exp(-1) ) ./ (1 - exp(-min_dist(too_close)/repulsion.d)) .* dir .* (1 - 2*is_inside(too_close));
% % if is_inside = 0, then we multiply by 1; if is_inside = 1, then we multiply by -1, meaning we actually have a force pulling the offending network node toward the
% % mesh, in the outward direction in hopes it will exit soon
%%

[ A_temp, A_force_recycle, A_torque_recycle, RHS, A_motor_torque0] = matrix_assembly_mex_wrapper(Mesh,Network, Repulsion, matrix_props,index_mapping,mesh_node_parameters,assembly_input, t);

[solution] = matrix_solve(assembly_input, matrix_props, A_temp, RHS);

der_fixed = solution( end-6 : end);  % since mesh was updated above, this is now in the fixed frame

% A_1 = A_1_matrix(y_swimmer(4:6)); % see comments in this function for detailed explanation of what the angles y(4:6) represent


% der_bodyframe(1:3) = A_1' * der_bodyframe(1:3);  % go from fixed frame to body frame using inverse or transpose of A_1, which was used to go from body to fixed
% der_local(4:6) = A_1' * der_local(4:6);
% now we've rotated U, Omega from fixed frame to body frame, so old algorithm based on body frame U, Omega applies?


 der = der_fixed; %initialization, will overwrite most values but will keep final value for omega of tail
 
% der(1:3) = A_1 * der_local(1:3);  % in the Ramia papers
% this gets cancelled by above A_1' * der_bodyframe(1:3)

% % see Maple sheet timestepping_eqs_omega
% % der(4:6) =   [sin(y(6))/sin(y(5))  ,  cos(y(6))/sin(y(5))  ,  0 ; ...
% %     cos(y(6))            ,  -sin(y(6))           ,  0 ; ...
% %     -cot(y(5))*sin(y(6)) ,  -cot(y(5))*cos(y(6)) ,  1]  * der_body(4:6);  % in Ramia papers but beware typo/inconsistencies there

% der(4:6) = [ 0               ,       sin(y_swimmer(6))/cos(y_swimmer(5))  ,       cos(y_swimmer(6))/cos(y_swimmer(5))  ; ...
%     0               ,       cos(y_swimmer(6))            ,       -sin(y_swimmer(6))           ; ...
%     1               ,       sin(y_swimmer(6))*tan(y_swimmer(5))  ,       cos(y_swimmer(6))*tan(y_swimmer(5))  ; ]         * der_local(4:6);

% this is the version for Omega in the fixed frame, see Goldstein for eq A.15xyz giving Omega as function of angles and d(angles)/dt and see Maple sheet for
% solving those for d(angles)/dt
phi = y_swimmer(4); theta = y_swimmer(5); psi = y_swimmer(6);
der(4:6) = [ sin(theta)*cos(phi)/cos(theta)               ,       sin(theta)*sin(phi)/cos(theta)  ,       1  ; ...
    -sin(phi)               ,       cos(phi)            ,       0          ; ...
    cos(phi)/cos(theta)              ,       sin(phi)/cos(theta)  ,      0  ; ]         * der_fixed(4:6);

% if abs(y(5)) <  1E-10 % 2nd angle is approaching 0 and der will blow up
%     der(4:6) = NaN;
% end

if abs(tan(y_swimmer(5))) >  1E10 % 2nd angle is approaching 0 and der will blow up
    %     der(4:6) = NaN;  %hopefully enough to signal a problem
    error('ODE blow up imminent');
end

% solution_fixedframe = [ reshape(  A_1 * reshape( solution(1:end-1), 3,[]) , [], 1); solution(end) ];

[u_network_nodes] = field_velocity(Mesh,Network, Repulsion, Network.nodes, solution, matrix_props,index_mapping,mesh_node_parameters,assembly_input);


% y = [l; u_network; U; Omega; omega]


dldt =  Network.E .* Network.l_0 ./ Network.eta .* ( r ./ Network.l - 1);


dydt = [ dldt; reshape(u_network_nodes', Network.n_nodes*3,1); der ];
