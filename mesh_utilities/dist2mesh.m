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
options = optimoptions('fmincon','OptimalityTolerance',1E-7,'ConstraintTolerance',1E-12,'StepTolerance',1E-9,...
    'MaxFunctionEvaluations',1E4,'MaxIterations',1E4,'Display','off');
xi_eta_edge_tol = 0.05;  %2E-3; % if we are within this dimensionless distance from any of the 3 sides of the triangle in xi, eta space, then we have to check adjacent element
% 1E-3 leads to rare failures to find the most correct normal vector, but
% will leave for now and instead relax tolerance on dotprod near end of
% code, shouldn't matter (?)

x = NaN(size(points,1),3);
xi_eta = NaN(size(points,1),2);
mesh_index = NaN(size(points,1),1);
element_index = NaN(size(points,1),1);
is_inside = false(size(points,1),1);
distance2mesh = NaN(size(points,1),1);

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

% interrogation_count = 0;
% interrogations = [];


for p = 1:size(points,1) % parfor is slower ?
    
    %     Dist2 = Inf(length(Mesh),1); % closest dist2 to each submesh, in case point is within danger zone of more than one
    solutions = NaN(0,5);
    %     is_closer = true(1,length(Mesh));
    for j = 1:length(Mesh)
        
        if is_close(p,j)
            
            %             [min_dist2_node, node_ind] = min(   sum(  (Mesh(j).nodes - points(p,:)).^2 , 2)  ); % index of closest node to x0
            
            [min_dist2_node, node_inds] = mink(   sum(  (Mesh(j).nodes - points(p,:)).^2 , 2) , 4 ); % indices of closest 4 nodes to x0
            
            
            if  min_dist2_node(1) > min_dists2(j,2)   % actually, this point is not that close to the nearest node, so change it to "not close" and continue
                %                 is_closer(j) = false;
                continue
            end
            
            %             interrogation_count = interrogation_count + 1; % how many times we've looked harder at distance to a point
            %             optcount = 0;
            % this is probably broken for geometries with coincident
            % submeshes that join at sharp corners - would have to do the
            % below over the combined set of elements instead of one
            % submesh a time ?
            
            % the method may also fail in cases of really complex meshes (e.g. donut-like curved rods) or really concave regions, where the closest element is actually distant
            % from the closest node in terms of connectivity (?)
            
            % not sure about above anymore now that code has been improved,
            % maybe no longer problems
            
            element_inds = []; % all elements we've interrogated for distance to current point
            
            node_info = { index_mapping.global_node2local_node{  index_mapping.local_node2global_node{j}(node_inds)  } };
            shared_elements = intersect( node_info{1}(:,3), node_info{2}(:,3) );
            switch length(shared_elements)
                case 1 % these nodes "span" across a single element, surely the closest point is within this element (?)
                    element_inds(end+1,1) = shared_elements;
                    
                case 2 % these 2 nodes are along an edge, can't be sure what element has closest point
                    
                    shared_elements_3 = intersect( shared_elements, node_info{3}(:,3) ); % are there any elements containing all 3 closest nodes?
                    
                    switch length(shared_elements_3)
                        
                        case 1 % narrowed it down to one element
                            element_inds(end+1,1) = shared_elements_3;
                        case 0 % failsafe, shouldn't happen?
                            element_inds(end+1,1) = shared_elements(1);
                        otherwise
                            shared_elements_4 = intersect( shared_elements_3, node_info{4}(:,3) ); % are there any elements containing all 4 closest nodes?
                            switch length(shared_elements_4)
                                case 0 % failsafe, shouldn't happen?
                                    element_inds(end+1,1) = shared_elements_3(1);
                                otherwise
                                    % if more than one, randomly pick the first shared element
                                    % containing all 4 nodes, hopefully we'll
                                    % find the correct one eventually
                                    element_inds(end+1,1) = shared_elements_4(1);
                            end
                    end
                    
                case 0
                    element_inds(end+1,1) = node_info{1}(1,3);
                    % if the closest 2 nodes don't share an element,
                    % take closest node and arbitrarily the first element
                    % it's part of and hope things converge
                otherwise
                    % this is super weird, worth investigating if it
                    % happens
                    error('closest 2 nodes share more than 2 elements, WTF');
                    %                    element_inds(end+1,1) = shared_elements(1);
            end
            
            xi_eta_temp = NaN(1,2);  dist2 = NaN(1,1);
            
            for c = 1:25 % if we aren't done after looking at 25 elements, something prolly wrong
                
                objfun = @(xe) dist(xe,points(p,:)',Mesh(j).nodes(Mesh(j).elements(element_inds(end),:), :),Mesh(j).shape_parameters(element_inds(end),:));
                
                flag = NaN;
                 rng(0);
                for g = 1:25 % if none of 25 random initial guesses don't work, something prolly wrong
                    if g == 1
                        guess = [1/3 1/3]'; % original choice seems to nearly always work....
                    else
                       
                        while true
                            guess = rand(2,1); % xi, eta between 0 and 1
                            if sum(guess) < 1 % within reference triangle
                                break
                            end
                        end
                    end
                    [xi_eta_temp(c,:), dist2(c,1), flag, output] = fmincon(objfun,guess,[1 1],1, [],[],[0 0]',[1 1]',[],options);  % linear constraint is that eta = 0 .. 1 - xi or 1*eta + 1*xi = 1
                    %                 output.message
                    if flag > 0
                        break
                    end
                end
                % exitflag % always seems to converge based on first-order optimality
                % (optimality tolerance)
                if flag <= 0
                    error('failure of fmincon');
                end
                % if output.bestfeasible.fval ~= dist2(c,1)
                %     pause
                % end
                %                 output.algorithm % seems to always use interior point algorithm
                %                 optcount = optcount + 1;
                % if we're actually close to one of the corner
                % vertices, we may choose the wrong adjacent element
                % below but we should eventually end up at the correct
                % element after additional iterations (?)
                if xi_eta_temp(c,1) < xi_eta_edge_tol % too close to left edge
                    midpt_node_ind = Mesh(j).elements(element_inds(end),6); % midpt node shared with adjacent element and no others
                elseif xi_eta_temp(c,2) < xi_eta_edge_tol % too close to bottom edge
                    midpt_node_ind =  Mesh(j).elements(element_inds(end),4); % midpt node shared with adjacent element and no others
                elseif  sum(xi_eta_temp(c,:)) > (1 - xi_eta_edge_tol) % too close to diagonal edge
                    midpt_node_ind = Mesh(j).elements(element_inds(end),5); % midpt node shared with adjacent element and no others
                else % not close to any edges, so this element presumably really contains the closest point on the mesh
                    break
                end
                
                
                
                node_info =  index_mapping.global_node2local_node{  index_mapping.local_node2global_node{j}(midpt_node_ind)  } ; %various index info for all (2) elements containing the midpt node
                node_info = node_info(node_info(:,1) == j,:); % only keep matches on current submesh
                element_ind = setdiff(node_info(:,3), element_inds(end)); % element index for "other" element containing this midpt node that we haven't just looked at
                % if the midpt node is actually on the edge of a 2D sheet,
                % element_ind will be empty since there is no adjacent element
                
                if isempty(element_ind) || ismember(element_ind, element_inds) % no adjacent element, or, the next element is one we've already looked at, so we're done (closest point on mesh is very close to an edge or vertex)
                    
                    break
                else
                    element_inds(end+1,1) = element_ind;
                end
                
            end % for cc = 1:25
            
            if c == 25
                error('max element iteration reached')
            end
            
            %             if c ~= 1
            %                 disp(['interrogated ',num2str(c),' elements']);
            %             end
            
            %             [min_dist2,ind] = min(dist2);  % element with the smallest distance to the point
            %             xi_eta_temp = xi_eta_temp(ind,:);
            %             element_ind = r(ind);
            
            
            solutions(end+1:end+length(element_inds),:) = [dist2, xi_eta_temp, repmat(j,length(element_inds),1), element_inds];
            
        end % is_close
        
    end % submeshes
    
    
    
    if ~isempty(solutions) % we've bothered computed distance from current point to at least one element
        
        [~,inds] = sort(solutions(:,1)); % element with closest distance is presumably correct one, but keep others around since if we are slightly wrong and near a sharp corner, the normal could be catestrophically incorrect
        distance2mesh(p) = sqrt(solutions(inds(1),1));
        mesh_index(p) = solutions(inds(1),4);
        element_index(p) = solutions(inds(1),5);
        xi_eta(p,:) = solutions(inds(1),2:3);
        
        for cc = 1:size(solutions,1)
            
            
            [x_temp,~,~,n] = T6interp(  Mesh(solutions(inds(cc),4)).nodes(    Mesh( solutions(inds(cc),4) ).elements(solutions(inds(cc),5),:),   :    ), ...
                solutions(inds(cc),2),     solutions(inds(cc),3) ,...
                Mesh( solutions(inds(cc),4) ).shape_parameters(solutions(inds(cc),5),:) );
            
            if cc == 1 % this was the best solution, output the position
                x(p,:) = x_temp;
                pvec = (points(p,:) - x(p,:));
            end
            
            
            
            dotprod = (pvec / sqrt(sum(pvec.^2))) * n;
            
            if abs(dotprod) > 0.5 % a good choice here depends on the sharpest corners of the mesh.  For 90 degree corners, we expect the worst case dotprod for the "correct" element's normal to be around 0.7 = cos(pi/4) for an outside point.
                % for sharper corners, the worst case dotprod for an outside point will decrease and
                % get close to zero (but still positive) as the angle of the
                % corner goes to zero.
                break % we have found a reasonable normal, don't need to look at other nearby elements
                
                
            end
            
            
            
        end % cc = 1:size(solutions,1)
        
        if abs(dotprod) < 0.3 %&&  distance2mesh(p) >
            %             error('could not find reasonable normal vector in any interrogated elements')
            %
            %             figure(834);  clf; [s,e] = plot_mesh(Mesh(1),3);
            %
            %             hold on
            %             nn = plot3(points(p,1),points(p,2),points(p,3),'ko','MarkerFaceColor','k');
            %             cp = plot3(x_temp(1),x_temp(2),x_temp(3),'ro','MarkerFaceColor','r');
            %             q1 = quiver3(x_temp(1),x_temp(2),x_temp(3),n(1),n(2),n(3),'k','LineWidth',1.5);
            %             pv = (pvec / sqrt(sum(pvec.^2)));  q2 = quiver3(x_temp(1),x_temp(2),x_temp(3),pv(1),pv(2),pv(3),'r','LineWidth',1.5);
            %             hold off
            
            % now that we're using a threshold for distance inside, just going to
            % assume that we are outside if we can't figure it out
            % expect this to happen often right after we put nodes back on the mesh,
            % since distance2mesh will be tiny and the vector from x to the node will
            % be nearly tangent to the mesh, and the dotprod nearly zero - just can't
            % be sure whether we're inside or outside but it doesn't really matter
            
            % leave is_inside(p) as false, the initialized value
            continue
        end
        
        
        
        if dotprod < 0
            is_inside(p) = true;
            %warning(["network node inside by " + num2str(distance2mesh(p))])
            
            %               is_inside(p) = 0;  % initialized as false
        end
        
    end
    
end


if any(is_inside)
    1;
end
% optcount