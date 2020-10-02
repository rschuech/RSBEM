function [Mesh] = rotateTail(Mesh, node_parameters,index_mapping,tail_ind, angle)
% rotates bacterial tail mesh around it's own (motor) axis by angle
% to allow for tails located and oriented anywhere, we first shift to origin, rotate, then shift back

     temp = Mesh(tail_ind).refpoints(:,1);
    Mesh(tail_ind) = shiftMesh(Mesh(tail_ind), -Mesh(tail_ind).refpoints(:,1));  % shift so that refpoint is at origin to do rotation around motor axis
    Mesh = rotateMesh(Mesh,node_parameters,index_mapping,tail_ind, angle  ,   Mesh(tail_ind).orientation(:,1) );
    Mesh(tail_ind) = shiftMesh(Mesh(tail_ind), temp);  % shift back to where it was
