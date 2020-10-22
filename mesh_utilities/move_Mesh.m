function [Mesh,mesh_node_parameters, rotmat] = move_Mesh(Mesh, mesh_node_parameters, index_mapping, movements)
%this simple wrapper is meant to take original mesh (located at origin and not translated
%or rotated yet) and to update it's position and orientation to the current
%values via translation and rotations around x,
%y, z axes

%movements = [    x y z                  Theta_x Theta_y Theta_z]
%            [    translation shifts      orientation angle changes]


%combined x,y,z rotation - do while mesh is still centered at the origin so
%rotations easily done around x, y, z axes
% http://en.wikipedia.org/wiki/Rotation_matrix    section on General rotations

[Mesh , mesh_node_parameters, rotmat] = rotateMesh(Mesh,mesh_node_parameters, index_mapping,1:length(Mesh),  movements(4:6));

%then, translation
Mesh = shiftMesh(Mesh, movements(1:3));  %refpoint is initial position, doesn't change (so only input orig Mesh to this function!)
%shiftMesh takes care of refpoints, centroid, and Centroid as well as
%verts

