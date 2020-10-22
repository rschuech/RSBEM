function [refpoint, motor_orientation] = get_rotational_references(Mesh, assembly_input)

% returns reference point (refpoint) used when calculating rotational velocity Omega and, for rotating flagellum cases, the axis around which the tail
% rotates (motor_orientation)

tail_ind = 0;
for i = 1:length(Mesh)
    if "Tail" == Mesh(i).name  % [Mesh.name] not allowed for code generation, hence this loop
        tail_ind = i; break;
    end
end

if assembly_input.rotating_flagellum
  
      motor_orientation = Mesh(tail_ind).orientation(:,1);  % updated to allow for tails attached at an angle - motor rotation axis should always be the centerline axis of the tail
    refpoint = Mesh(tail_ind).refpoints(:,1); %make sure refpoint is somewhere along the axis of motor rotation for free swimming problems involving a motor torque condition
else
    motor_orientation = NaN(3,1);  %appease Coder
    refpoint = Mesh(1).refpoints(:,1); % doesn't really matter where it is, so use body or whatever first submesh is
end
