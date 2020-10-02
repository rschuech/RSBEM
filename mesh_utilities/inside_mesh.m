function [is_inside] = inside_mesh(points, Mesh)
% check whether each row of points is inside flat triangular version of
% Mesh (won't be exact for curved triangular meshes)

is_inside = inpolyhedron(Mesh.elems(:,1:3),Mesh.verts,points);  %way faster than my shite method
% 
% is_inside = NaN(size(points,1),1);
% 
% parfor i = 1:size(points,1);
%     pos = [];
%   %  disp(['inside mesh doneness = ',num2str(i/size(points,1))]);
%     point = points(i,:);
%     
%     n_intersections = 0;
%     while n_intersections == 0
%         rand_dir = 2 * rand(1,3) - 1;  %random [i j k] vector, each entry between -1 and 1
%         rand_dir = rand_dir / sqrt(sum(rand_dir.^2));  %normalize to make it a unit vector
%         LINE = [point rand_dir];
%         
%         [~,pos] = intersectLineMesh3d(LINE, Mesh.verts, Mesh.elems(:,1:3)); %only take flat corner verts if T6 Mesh was input
%         n_intersections = size(pos,1);
%     end
%     
%     left_intersections = sum(pos < 0);
%     right_intersections = sum(pos > 0);
%     
%     if rem(left_intersections,2) == 0 && rem(right_intersections,2) == 0 %even number intersections for both rays
%         is_inside(i) = false;
%     else
%         is_inside(i) = true;
%     end
%     
% end
% 
% if any(isnan(is_inside(:)))
%     disp('problemo')
%     pause
% end
