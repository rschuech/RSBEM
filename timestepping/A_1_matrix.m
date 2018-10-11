function [A_1] = A_1_matrix(angles)




% eq 15 in Numerical model for locomotion of spirilla or eq 4.47 in
% Classical Mechanics textbook (Goldstein)
% A_1 = [cos(angles(3))*cos(angles(1)) - cos(angles(2))*sin(angles(1))*sin(angles(3))  ,  -sin(angles(3))*cos(angles(1)) - cos(angles(2))*sin(angles(1))*cos(angles(3))  ,  sin(angles(2))*sin(angles(1)); ...
%     cos(angles(3))*sin(angles(1)) + cos(angles(2))*cos(angles(1))*sin(angles(3)) ,  -sin(angles(3))*sin(angles(1)) + cos(angles(2))*cos(angles(1))*cos(angles(3))  ,  -sin(angles(2))*cos(angles(1)); ...
%     sin(angles(2))*sin(angles(3))                                 ,  sin(angles(2))*cos(angles(3))                                   ,  cos(angles(2))            ];



% Classical Mechanics textbook (Goldstein)  XYZ convention in appendix
% angles = Goldstein's [phi  theta  psi]
% 2 equivalent interpretations of the rotations: (see wikipedia discussion
% of extrinsic, intrinsic rotations)
% intrinsic:
% first rotation of phi around z, next rotation of theta around current y, last rotation of psi around current x
% extrinsic:
% first rotation of psi around fixed x, next rotation of theta around fixed y, last rotation of phi around fixed z

% for intrinsic, rotations are done left to right (Z Y' X'') 
% for extrinsic, rotations are done right to left (Z Y X)
% results are the same

% so, extrinsic interpreation is fine, but remember that psi = x, theta = y, phi = z

% note typo in original textbook eq, corrected on errata website and here!

A_1 = [cos(angles(2))*cos(angles(1))  ,  sin(angles(3))*sin(angles(2))*cos(angles(1)) - cos(angles(3))*sin(angles(1))  ,  cos(angles(3))*sin(angles(2))*cos(angles(1)) + sin(angles(3))*sin(angles(1)); ...
       cos(angles(2))*sin(angles(1))  ,  sin(angles(3))*sin(angles(2))*sin(angles(1)) + cos(angles(3))*cos(angles(1))  ,  cos(angles(3))*sin(angles(2))*sin(angles(1)) - sin(angles(3))*cos(angles(1)); ...
       -sin(angles(2))        ,  cos(angles(2))*sin(angles(3))                         ,  cos(angles(2))*cos(angles(3))                       ;  ];