function [forces] = integrate_traction(A_force,A_torque,f)


%for i_mesh = 1:length(Mesh)
   % n_unq_verts = diff(Mesh(i_mesh).indices.glob.unq_bounds.vert) + 1;

   
   for i_mesh = 1:size(A_force,3)
       integrals = [A_force(:,:,i_mesh); A_torque(:,:,i_mesh)]  *  f;  %cleverly, each page of A_force and A_torque (going with each submesh) has zeros for all the verts not part of the current submesh, so they can each be multiplied by full f
   
   
forces(i_mesh).drag = -integrals(1:3);  %A_force_torque*f is force of surface on water; want force of water on surface
forces(i_mesh).torque = -integrals(4:6);

   end

% A_force_torque = [sum(A_force,3); sum(A_torque,3)]; %combine integrals over each submesh since we want force and torque over entire object
% integrals = A_force_torque * f;
% 
% 
