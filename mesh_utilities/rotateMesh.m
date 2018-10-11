function [Mesh, rotmat] = rotateMesh(Mesh,angles,varargin)
%rotates mesh in various ways:

% default is to do intrinsic rotations around z, new y, and new x axes as
% in Goldstein textbook appendix  (or equivalently, extrinsic rotations
% around fixed x, fixed y, fixed z)
% in both cases, angles = [Z Y X]
% varargin = {}

% or
% around x, y, z axes through origin by angles (extrinsic rotations)
% varargin = {[ax1 ax2 ax3} = order of rotations
% or
% around arbitrary vector by some angle            varargin = {vector}
% around z, y, x intermediate axes by angles (really just a scalar) (intrinsic rotations)

 %combined x,y,z rotation - do while mesh is still centered at the origin so
    %rotations easily done around x, y, z axes
    % http://en.wikipedia.org/wiki/Rotation_matrix    section on General rotations
    if nargin == 2 % extrinsic rotations in order z, y, x
       
       rotmat = A_1_matrix(angles);
    elseif length(angles) == 1 % have input a single angle and rotation vector to rotate around arbitrary vector
        rotmat = rotate_arbitrary_vector( varargin{1}, angles);
        %otherwise, have optionally input just a list of order indices
        %for extrinsic rotations around origin
    elseif isequal([varargin{1}(:)], [1 2 3]') % x y z
        rotmat = rotation_matrix('z',angles(3)) * rotation_matrix('y',angles(2)) * rotation_matrix('x',angles(1));   
    elseif isequal([varargin{1}(:)], [2 1 3]') % y x z
        rotmat = rotation_matrix('z',angles(3)) * rotation_matrix('x',angles(1)) * rotation_matrix('y',angles(2));
    elseif isequal([varargin{1}(:)], [1 3 2]') % x z y
        rotmat = rotation_matrix('y',angles(2)) * rotation_matrix('z',angles(3)) * rotation_matrix('x',angles(1));
    elseif isequal([varargin{1}(:)], [3 1 2]') % z x y
        rotmat = rotation_matrix('y',angles(2)) * rotation_matrix('x',angles(1)) * rotation_matrix('z',angles(3));
    elseif isequal([varargin{1}(:)], [3 2 1]') % z y x
        rotmat = rotation_matrix('x',angles(1)) * rotation_matrix('y',angles(2)) * rotation_matrix('z',angles(3));
    elseif isequal([varargin{1}(:)], [2 3 1]') % y z x
        rotmat = rotation_matrix('x',angles(1)) * rotation_matrix('z',angles(3)) * rotation_matrix('y',angles(2));
    else
        rotmat = NaN(3,3); %something went wrong....
    end
    
    
for i = 1:length(Mesh)
    
   
        
    Mesh(i).verts = (rotmat * Mesh(i).verts')';
    %also update vectors and other coordinate variables
    Mesh(i).orientation = rotmat * Mesh(i).orientation; %orientation vectors only need to be rotated, not shifted
    Mesh(i).orientation = Mesh(i).orientation ./ repmat(  sqrt(sum(Mesh(i).orientation.^2)),3,1);
    
    Mesh(i).refpoints = rotmat * Mesh(i).refpoints;
    Mesh(i).centroid = rotmat * Mesh(i).centroid;
    Mesh(i).Centroid = rotmat * Mesh(i).Centroid;
    
    
    
end