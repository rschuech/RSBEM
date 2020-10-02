function [rotmat] = rotation_matrix(axis,theta)

switch axis
    case 'x'
        rotmat = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    case 'y'
        rotmat = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    case 'z'
        rotmat = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
end
