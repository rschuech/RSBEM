function [is_intersected] = submesh_intersections(Mesh,node_parameters,index_mapping, dist_tol ,n_angles, plot_intersections, nthreads,varargin)
%checks for collisions between body and tail during tail rotation, which would mean an
%invalid geometry.  n_angles evenly spaced phase angles are checked.
% varargin is the mesh refinement levels
if isempty(varargin)
    nrefines = [3 2];
else
    nrefines = varargin{1};
end
% disp('Checking for submesh intersections');


%dist_tol =  0.031018 / 2; %half the tail radius
if ~isempty(n_angles)   % this is prolly bacteria mesh, we need to check many tail phase angles by rotating tail mesh here
    angles = linspace(0,2*pi,n_angles); %phase angles to test
    is_intersected = false(1,length(angles));
    
    if nthreads > 1
        init_parpool(nthreads, length(angles));
    end
    
    % subsample meshes to decrease liklihood of missing an actual or near intersection
    Surface_refined = refine_vis_surface(Mesh, nrefines);  %refine body more than tail
    
    
%     parfor (i = 1:length(angles), nthreads)  %at least for 10 test angles, parfor seems about 6 times faster than for
        for i = 1:length(angles)
        angle = angles(i);
%         Mesh_temp = struct;
        Mesh_temp = Surface_refined(1);
%         [Mesh_temp(2)] = rotateMesh(Surface_refined(2), [   0  0  angle]' );  %rotate tail around x axis
         [Mesh_temp] = rotateTail(Surface_refined, node_parameters,index_mapping,2, angle);
        %     idx = rangesearch(Mesh_temp(1).nodes,Mesh_temp(2).nodes,dist_tol); %finds all the X points that are within distance r of the Y points. Rows of X and Y correspond to observations, and columns correspond to variables.
        %     intersections = [];
        %     for id = 1:length(idx)
        %         if ~isempty(idx{id})
        %             intersections(end + 1 : end + length(idx{id}),:) = Mesh_temp(1).nodes(idx{id},:);
        %
        %         end
        %     end
        
        % this choice of tail and then body input order seems slightly faster than above
        % choice of body and then tail
        idx = rangesearch(Mesh_temp(2).nodes,Mesh_temp(1).nodes,dist_tol); %finds all the X points that are within distance r of the Y points. Rows of X and Y correspond to observations, and columns correspond to variables.
        intersections = [];
        for id = 1:length(idx)
            if ~isempty(idx{id})
                intersections(end + 1 : end + length(idx{id}),:) = Mesh_temp(2).nodes(idx{id},:);
                
            end
        end
        
        % intersections = mesh2mesh(Mesh(1).elems(:,1:3),Mesh(1).nodes,MeshR.elems(:,1:3),MeshR.nodes);
        
        if plot_intersections && ~isempty(intersections)
            
            
            figure(635);  cla;  s = plot_mesh( rotateTail(Mesh,node_parameters,index_mapping,2, angle), [1 1]); hold on;  set(s,'facealpha',1);  light;
            if ~isempty(intersections)
                ph = plot3(intersections(:,1),intersections(:,2),intersections(:,3),'o','markerfacecolor','b','markersize',8);
            end
            title(['phase angle = ',num2str(angle)]);
            size(intersections,1)
            drawnow
            if ~isempty(intersections)
                pause
            end
        end
        
        if ~isempty(intersections)
            is_intersected(i) = true;
            
            %break %don't bother checking any other phase angles if we're already screwed
            
            %ordinarily not including the break would be pretty crappy (break
            %not allowed in parfor) but if length(angles) < num CPUs, it
            %shouldn't really matter since all angles are tested simultaneously
            %and a parfor iteration can't be stopped mid-computation
            %(i.e. 10 angle tests should take the same amount of time as 1
            %angle test since they're all done at once)
        end
        
    end
    
    
    is_intersected = any(is_intersected);  %if it ever intersected, it's no good
    
else  %prolly dinoflagellate - only need to check for intersections between body, tail in the exact mesh that was input
    
    % subsample meshes to decrease liklihood of missing an actual or near intersection
    
%     if length(Mesh) == 3 %no wingtip
%         nrefines = [3 0 1];
%     elseif length(Mesh) == 4
%         nrefines = [3 0 1 0];
%     end
    
    nrefines = 1;
    
    Surface_refined = refine_vis_surface(Mesh, nrefines);  %refine body more than tail; don't need to refine transverse at all since not checking it
    
   % Mesh_temp = struct;
    tail_ind = find(strcmp({Mesh.name} , 'Tail'));
    other_inds = setdiff(1:length(Mesh),tail_ind);
    
%     Mesh_temp = Surface_refined(other_inds);  %should always be body, transverse, hairs
%     Mesh_temp(2) = Surface_refined(3); %should always be tail....
    tail_nodes = Mesh(tail_ind).nodes;
    other_nodes = vertcat(Mesh(other_inds).nodes);
    
    idx = rangesearch(tail_nodes,other_nodes,dist_tol); %finds all the X points that are within distance r of the Y points. Rows of X and Y correspond to observations, and columns correspond to variables.
    
    intersections = [];
    for id = 1:length(idx)
        if ~isempty(idx{id})
            intersections(end + 1 : end + length(idx{id}),:) = tail_nodes(idx{id},:);
            
        end
    end
    
    if plot_intersections && ~isempty(intersections)
        figure(635);  cla;  s = plot_mesh( Mesh, [2]); hold on;  set(s,'facealpha',0.7);  light;
        
        ph = plot3(intersections(:,1),intersections(:,2),intersections(:,3),'o','markerfacecolor','b','markersize',8);
        
        size(intersections,1)
        drawnow
        if ~isempty(intersections)
            pause
        end
    end
    
    
    
    if ~isempty(intersections)
        is_intersected = true;
    else
        is_intersected = false;
    end
    
    
    
end
