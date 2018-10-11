function rotmat = rotate_arbitrary_vector( vec, theta)

% rotation matrix for rotating around an arbitrary vector starting at
% origin by angle

% vec better be a unit vector
vec = vec / sqrt(sum(vec.^2));
% http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/

rotmat = [ vec(1)^2 + (1 - vec(1)^2) * cos(theta)   ,    vec(1)*vec(2)*(1-cos(theta)) - vec(3)*sin(theta)   ,    vec(1)*vec(3)*(1-cos(theta)) + vec(2)*sin(theta) ; ...
           vec(1)*vec(2)*(1-cos(theta)) + vec(3)*sin(theta)   ,   vec(2)^2 + (1 - vec(2)^2)*cos(theta)    ,    vec(2)*vec(3)*(1 - cos(theta)) - vec(1)*sin(theta) ; ...
           vec(1)*vec(3)*(1 - cos(theta)) - vec(2)*sin(theta)  ,  vec(2)*vec(3)*(1 - cos(theta)) + vec(1)*sin(theta)  ,   vec(3)^2 + (1 - vec(3)^2)*cos(theta) ;  ];
       
       