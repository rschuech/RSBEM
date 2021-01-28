function [distance2mesh, x, xi_eta, mesh_index, element_index, is_inside] = dist2mesh(points, Mesh, index_mapping, min_dists2)
% points is N x 3 to check closest distance to all submeshes, when distance is within some danger zone
% min_dists2 is n_submeshes x 2, first column is min distance^2 away from body centroid or tail axis line segment to be worth caring about
% 2nd column is min distance^2 away from closest mesh node to be worth caring about
% if point passes both those tests, we compute the exact distance to the closest point on the closest mesh element

% x is N x 3 closest point on mesh
% xi_eta N x 2 are xi and eta for x
% element_index N x 1 is which local element index we're talking about
% is_inside N x 1 is whether the point is inside the mesh or not
% all outputs are NaN if we didn't bother computing exact distance since the point was not within the danger zone
options = optimoptions('fmincon','OptimalityTolerance',1E-9,'ConstraintTolerance',1E-12,'StepTolerance',1E-9,'Display','notify');

x = NaN(size(points,1),3);
xi_eta = NaN(size(points,1),2);
mesh_index = zeros(size(points,1),1);
element_index = NaN(size(points,1),1);
is_inside = false(size(points,1),1);
distance2mesh = zeros(size(points,1),1);

is_close = false(size(points,1),length(Mesh));

for i = 1:length(Mesh)
    
    switch Mesh(i).name
        case 'Body'
            
            is_close(:,i) = sum(  (Mesh(i).refpoints(:,2)' - points).^2 , 2) < min_dists2(i,1);
            
            
        case 'Tail'
            % http://geomalgorithms.com/a02-_lines.html
            % assume that tail axis line segment starts at first refpoint (near body) and ends at 2nd refpoint, which should correspond to distal end of tail
            w = points - Mesh(i).refpoints(:,1)';
            v_L = Mesh(i).refpoints(:,2) - Mesh(i).refpoints(:,1);
            b = w * v_L / (v_L' * v_L); % value of line parameter corresponding to closest point on line to each points coord
            before_start = b < 0;
            after_end = b > 1;
            in_between = ~before_start & ~after_end;
            is_close(before_start,i) = sum(  (Mesh(i).refpoints(:,1)' - points(before_start,:)).^2 , 2) < min_dists2(i,1);
            is_close(after_end,i) = sum(  (Mesh(i).refpoints(:,2)' - points(after_end,:)).^2 , 2) < min_dists2(i,1);
            is_close(in_between,i) = sum( (  w(in_between,:) - b(in_between).*repmat(v_L',sum(in_between),1)  ).^2 , 2) < min_dists2(i,1);
    end
end


for p = 1:size(points,1) % parfor is slower ?
    
    %     Dist2 = Inf(length(Mesh),1); % closest dist2 to each submesh, in case point is within danger zone of more than one
    solutions = NaN(length(Mesh),4);
    is_closer = true(1,length(Mesh));
    for j = 1:length(Mesh)
        
        if is_close(p,j)
          
%             [min_dist2_node, node_ind] = min(   sum(  (Mesh(j).nodes - points(p,:)).^2 , 2)  ); % index of closest node to x0
         
            [min_dist2_node, node_ind] = mink(   sum(  (Mesh(j).nodes - points(p,:)).^2 , 2) , 2 ); % indices of closest 2 nodes to x0
                  
            
            if min_dist2_node > min_dists2(j,2)  % actually, this point is not that close to the nearest node, so change it to "not close" and continue
                is_closer(j) = false;
                continue
            end
            
            
            [r,~] = find(node_ind == Mesh(j).elements);  % indices of all elements containing this node
            % make decent guess that closest element to x0 is one of the elements containing the closest node to x0
            % should only be violated in cases of really complex meshes (e.g. donut-like curved rods) or really concave regions, where the closest element is actually distant
            % from the closest node in terms of connectivity (?)
            
            
            xi_eta_temp = NaN(length(r),2);  dist2 = NaN(length(r),1);
            
%           tic
            for c = 1:length(r)
                
                
                objfun = @(xe) dist(xe,points(p,:)',Mesh(j).nodes(Mesh(j).elements(r(c),:), :),Mesh(j).shape_parameters(r(c),:));
                
                [xi_eta_temp(c,:), dist2(c,1)] = fmincon(objfun,[1/3 1/3]',[1 1],1, [],[],[0 0]',[1 1]',[],...  % linear constraint is that eta = 0 .. 1 - xi or 1*eta + 1*xi = 1
                    options);
            end
            
            %  right here, don't yet throw out all solutions but the
            %  closest.  Check if any solutions have nearly the same
            %  distance.  If so, keep all those around.  Do following calcs
            %  for all of them.  The x and pvec will be basically the same for all,
            %  but the normal n could be very different if the closest
            %  point is near a corner.  Compute dotprod for all cases.
            %  Find the one with the largest abs magnitude, closest to -1
            %  or +1.  Others may be close to 0, those are not reliable.
            %  The one closest to -1 or +1 should be correct and tell us if
            %  it's inside or outside.
            
            % also find way to only calc over a few elements by finding k
            % nearest mesh nodes to the network node?
%          toc
         
            [min_dist2,ind] = min(dist2);  % element with the smallest distance to the point
            xi_eta_temp = xi_eta_temp(ind,:);
            element_ind = r(ind);
            
            
            solutions(j,:) = [min_dist2, xi_eta_temp, element_ind];
            
        end % is_close
        
    end % submeshes
    
    if any(  is_close(p,:) & is_closer    )
        [~,ind] = min(solutions(:,1));
        distance2mesh(p) = sqrt(solutions(ind,1));
        xi_eta(p,:) = solutions(ind,2:3);
        mesh_index(p) = ind;
        element_index(p) = solutions(ind,4);
        [x(p,:),~,~,n] = T6interp(Mesh(ind).nodes(Mesh(ind).elements(element_index(p),:), :), solutions(ind,2),solutions(ind,3)     ,  Mesh(ind).shape_parameters(element_index(p),:) );
        
        pvec = (points(p,:) - x(p,:));
        dotprod = (pvec / sqrt(sum(pvec.^2))) * n;
        
%         if abs( abs(dotprod) - 1 ) > 1E-3
%            
%             warning('Check min distance calc');  % this test is expected to fail if the point is near a mesh corner - if the closest point on the mesh is along an
%             % element edge, we should really be getting the average normal
%             % there, not just taking the normal from one of the adjacent
%             % elements - of course that will be wrong.  This code seems to
%             % be working otherwise, so not doing this check anymore.  The
%             % direction of the repulsive force is computed using the two
%             % points, not the mesh normal, so there's no problem even near
%             % corners. The is_inside test below should also still work, the
%             % dotprod just won't always be nearly 1 or -1 but it should
%             % have the correct sign (unless there's a *really* sharp convex
%             corner or the point is *really* in the middle of two
%             elements?)
%         end
        
        if dotprod < 0
            is_inside(p) = true;
            warning(["network node inside by " + num2str(distance2mesh(p))])
        else
%               is_inside(p) = 0;  % initialized as false
        end
        
    end
    
end