


% if Metadata(1).geom.AR2 == 0 %straight rod, origin is in center of body
%     body_rear_pole = [0 - Metadata(1).geom.totlength / 2  0   0];
% else
%     body_rear_pole = [0 0 0]' - [Metadata(1).geom.radius 0 0]';
% end
% 
% tail_front_pole = Mesh(2).refpoints + [Metadata(2).geom.pipeRadius 0 0]';
% 
% midpoint = (body_rear_pole(:) + tail_front_pole(:)) / 2;  %more or less middle of hook
% 
% Mesh(2).refpoints(:,1) = midpoint;  %this is what integrate_moments uses to calculate torque

% tiny difference between original tail refpoint (center of end sphere) and
% middle of "hook" won't be significant - annoying to calculate
% latter for pole2pole configuration


clear force torque torque_ratio
torques = NaN(3,length(Solutions.f));

for fi = 1:length(Solutions.f)
    
    Mesh_temp = Mesh;
    
    y = [Mesh(1).refpoints(:,1); [0 0 0 Solutions.phase(fi)]' ];
    
    
    %first rotate only tail to handle tail phase angle relative to body
%     Mesh_temp(2) = rotateMesh(Mesh_temp(2),[  0 0 y(7)]');  %rotate tail by y(7) around x-axis
    [Mesh_temp(2) ] = rotateTail(Mesh_temp(2) , y(7));
    %then move all submeshes according to translation and rotation of body +
    %tail
%     Mesh_temp = move_Mesh(Mesh_temp,[y(1:3) - Mesh_temp(1).refpoints(:,1); y(4:6)]); % shift ends up being zero
    %yprime keeps track of body refpoint by convention, so define translation shift using original body refpoint
    
    
    
    [VL_unity, A_force] = compute_force_integral(Mesh_temp, matrix_props, input);
    [VL_moments, A_torque] = compute_torque_integral(Mesh_temp, matrix_props, input);
    
    
    f = Solutions.f{fi};
    
    clear integrals forces
    for i_mesh = 1:size(A_force,3)
        integrals = [A_force(:,:,i_mesh); A_torque(:,:,i_mesh)]  *  f;  %cleverly, each page of A_force and A_torque (going with each submesh) has zeros for all the verts not part of the current submesh, so they can each be multiplied by full f
        
        
        forces(i_mesh).drag = -integrals(1:3);  %A_force_torque*f is force of surface on water; want force of water on surface
        forces(i_mesh).torque = -integrals(4:6);
        
    end
    
    if any(abs(sum(horzcat(forces.drag),2)) > eps*100) || any(abs(sum(horzcat(forces.torque),2)) > eps*100)
        stopafra
    end
    
    % take component of drag force in direction of tail centerline for
    % buckling force and resultant torque in plane normal to tail v zzzz
    % centerline for bending torque (not worrying that the direction
    % probably changes with phase angle)
    
    % only works for original tail orientation
%     force(fi) = abs(forces(1).drag(1));  
%     torque(fi) = sqrt( forces(1).torque(2)^2 + forces(1).torque(3)^2 );
    
    rotmat = Mesh(2).orientation'; % left-multiplying by new basis vectors (in row form!) yields vector in new coord system  (equivalent to taking inverse of new basis vectors in column form)
    
    temp = rotmat * forces(1).drag; % all drag components
    force(fi) = abs(temp(1));  % component of force along tail axis (the other two are shear forces on hook I guess?)
    temp = rotmat * forces(1).torque; % all torque components
    torque(fi) = sqrt(temp(2)^2 + temp(3)^2);
    
%     torque_max(fi) = max(abs(temp(2:3)));
%     torque_total(fi) = sqrt(sum(temp.^2));
    
    torque_ratio(fi) = abs( temp(1) / torque(fi) );
    
    torques(:,fi) = temp;
end

% torques



if any( abs(torques(1,:) - input.tail.motor_torque) > eps )
    stopa
end

    % rescale force and torque, normalized to constant power for all bugs,
    % as was done for swimming speed
force = sqrt( force.^2 .* input.constants.power  ./ Results(r).Avg_Power);
torque = sqrt( torque.^2 .* input.constants.power  ./ Results(r).Avg_Power);
% don't need to do any rescaling of a ratio since it won't change
   
 [interpolant] = trig_interp_fit(Solutions.phase,force);
  fun = @(x) trig_interp_eval(interpolant, x); 
 buckling.force.mean = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
buckling.force.max = max(trig_interp_eval(interpolant, phases));  %dumb brute force method but should be accurate enough
 
 [interpolant] = trig_interp_fit(Solutions.phase,torque);
  fun = @(x) trig_interp_eval(interpolant, x); 
 buckling.torque.mean = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
buckling.torque.max = max(trig_interp_eval(interpolant, phases));  %dumb brute force method but should be accurate enough

 [interpolant] = trig_interp_fit(Solutions.phase,torque_ratio);
  fun = @(x) trig_interp_eval(interpolant, x); 
 buckling.torque_ratio.mean = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
buckling.torque_ratio.max = max(trig_interp_eval(interpolant, phases));  %dumb brute force method but should be accurate enough

    
    