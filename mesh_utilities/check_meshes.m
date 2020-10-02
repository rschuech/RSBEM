%%

nrefines = 1;
volume_tol = 0.005;  % <= 0.5% error allowable
%area_outlier_tol = 10;
max_angle = 90*pi/180 ; %max angle between adjacent element normals, located at centroids of the elements

skip_succeeded = false; % skip meshes with Metadata files already marked as succeeded

plot_bads = true;
plot_goods = true;

% unrefine crack
% [1.125 0.35; 1.25 0.375; 1.375 0.4; 1.375 0.425; 1.5 0.45; 1.625 0.475; 1.75 0.525; 1.875 0.575; 2 0.6; 2.25 0.65; 2.25 0.675; 
% unrefine crack A LOT
% [1.5 0.475; 1.625 0.5; 2 0.625; 2.25 0.7
% refine crack
% [1.375 0.3;

input.performance.nthreads = 10;

input.accuracy.integration_tol.area.abstol = 1E-9 ;
input.accuracy.integration_tol.area.reltol = 0.1;
input.accuracy.integration_tol.area.maxevals = Inf;

input.accuracy.integration_tol.centroid.abstol = 1E-9;
input.accuracy.integration_tol.centroid.reltol = 0.1;
input.accuracy.integration_tol.centroid.maxevals = Inf;

input.accuracy.integration_tol.volume.abstol = 1E-9;
input.accuracy.integration_tol.volume.reltol = 0.1;
input.accuracy.integration_tol.volume.maxevals = Inf;

input.performance.verbose = false;




if strcmp(geom(1).shape,'curved_rod') || strcmp(geom(1).shape,'tail')  %need to check each curved rod mesh for self intersection
    
    succeeded = false(1,length(geom));
    reason = NaN(1,length(geom));
    
    
    for c = 1:length(geom)  
        
        if c < 0
            continue
        end
%         c
        %                 if succeeded(c)
        %                     continue
        %                 end
        
      %  c/length(geom)
        
        path_ind = c;  %must be c for things to work since path variable contains full path including name
        
        meshname = paths(path_ind).mesh_file;
        [Mesh, Metadata] = load_mesh(paths(path_ind).mesh_file);  %all paths should be the same for a given shape
        
        if isempty(Metadata)
            disp('empty Metadata')
%             pause
        else
            if skip_succeeded && isfield(Metadata.mesh, 'meshing_succeeded') && Metadata.mesh.meshing_succeeded
                succeeded(c) = true;
                continue
            end
        end
        
        
        
        
        if isempty(Mesh)
            if input.performance.verbose
            disp([paths(path_ind).mesh_file,'      never created']);
            end
%              pause
            Metadata.mesh.meshing_succeeded = false;
            save(paths(c).metafile,'Metadata');
            succeeded(c) = false;
            reason(c) = -1;
%              pause
            
            continue
            
        end
        
        
        %Metadata should always exist, maybe not Mesh
        [Mesh] = global_inds(Mesh);
        [Mesh] = renumber_Mesh(Mesh);
        %         clear temp_input
        temp_input.performance.nthreads = input.performance.nthreads;
        temp_input.accuracy.integration_tol.area = input.accuracy.integration_tol.area;
        temp_input.accuracy.integration_tol.centroid = input.accuracy.integration_tol.centroid;
        temp_input.accuracy.integration_tol.volume = input.accuracy.integration_tol.volume;
        [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
        
        edges = get_edges(Mesh);  %used for both topology check as well as adjacent normals check
        
        metric = size(unique(Mesh.elems(:,1:3)),1) - size(edges.verts,1) + Mesh.n_elem;  % # corner verts - # edges + # elements     = 2 for sphere, 0 for torus
        
        if metric ~= 2 %don't have a topological sphere, something is wrong, probably self-intersection
            if metric == 0
                disp([paths(path_ind).mesh_file,'    is a topological torus']);
            else
                disp([paths(path_ind).mesh_file,'    is topologically wierd']);
            end
            
            if plot_bads
                figure(35)
                cla
                [s,e] = plot_mesh(Mesh,nrefines);
                set(s,'facealpha',1);
                title(paths(path_ind).mesh_file,'interpreter','none');
                drawnow
                %   disp('Mark mesh as bad and continue?');
%                 pause
            end
            
            
            %update Metadata file with self intersection flag
            
            if metric == 0
                Metadata.mesh.topology = 'torus';
                reason(c) = 1;
            else
                Metadata.mesh.topology = 'other';
                reason(c) = 2;
            end
            
            succeeded(c) = false;
            Metadata.mesh.meshing_succeeded = false;
            save(paths(path_ind).metafile,'Metadata');
            continue
            
        else %mesh appears OK
            Metadata.mesh.topology = 'sphere';
            %wait, mesh still may be broken in other ways
            %             Metadata.mesh.meshing_succeeded = true;
            %             succeeded(c) = true;
            %             reason(c) = 0;
        end
        
        % not out of the woods yet; may still have the wrong volume
        if ~strcmp(geom(1).shape,'tail') && (  Mesh.Volume < (geom(c).V*(1-volume_tol)) || Mesh.Volume > (geom(c).V*(1+volume_tol))  ) %if volume is either too small or too large
            if input.performance.verbose
            disp([paths(path_ind).mesh_file,'    has volume problem, error = ',num2str(100*(geom(c).V-Mesh.Volume)/geom(c).V),' %']);
            end
            if plot_bads
                figure(35)
                cla
                [s,e] = plot_mesh(Mesh,nrefines);
                set(s,'facealpha',1);
                title({paths(path_ind).mesh_file,['has volume problem, error = ',num2str(100*(geom(c).V-Mesh.Volume)/geom(c).V),' %']},'interpreter','none');
                drawnow
%                  pause
            end
            Metadata.mesh.meshing_succeeded = false;
            Metadata.mesh.wrong_volume = true;
            succeeded(c) = false;
            reason(c) = 3;
            save(paths(path_ind).metafile,'Metadata');
            continue
        end
        
        
        %       temp = abs((max(Mesh.area) -  mean(Mesh.area)) / std(Mesh.area));  %sunken elements hopefully have much larger area than mean, relative to stdev
        %          if temp > area_outlier_tol
        %             disp([paths(path_ind).mesh_file,'    has element outlier problem, # stdevs = ',num2str(temp)]);
        %               figure(34)
        %             cla
        %             [s,e] = plot_mesh(Mesh);
        %             set(s,'facealpha',0.4);
        %             title({paths(path_ind).mesh_file,['has element outlier problem, # stdevs = ',num2str(temp)]},'interpreter','none');
        %             drawnow
        %           %  pause
        %             Metadata.mesh.meshing_succeeded = false;
        %             Metadata.mesh.outliers = true;
        %             succeeded(c) = false;
        %             reason(c) = 3;
        %          end
        
        
        %check centroid normals of neighboring elements for jumps in direction
        normals = NaN(Mesh.n_elem,3);
        for i_elem = 1:Mesh.n_elem
            [~, ~, ~, normals(i_elem,:)] = T6interp(Mesh.verts(Mesh.elems(i_elem,:),:),0.5,0.5,Mesh.elem_params(:,i_elem));
        end
        
        angle = NaN(size(edges.verts,1),1);
        for i_edge = 1:size(edges.verts,1) %just need to loop over each edge, since each edge has two neighboring elements
            neighbor_normals = normals( edges.elems{i_edge}, :);
            angle(i_edge) = acos(dot(neighbor_normals(1,:),neighbor_normals(2,:)));
        end
        
        if max(angle) >= max_angle
            if input.performance.verbose
            disp([paths(path_ind).mesh_file,'    has normals angle problem']);
            end
            if plot_bads
                figure(35)
                cla
                [s,e] = plot_mesh(Mesh,nrefines);
                set(s,'facealpha',0.6);
                title({paths(path_ind).mesh_file,['has normals angle problem']},'interpreter','none');
                drawnow
                             pause
            end
            Metadata.mesh.meshing_succeeded = false;
            succeeded(c) = false;
            reason(c) = 4;
            save(paths(path_ind).metafile,'Metadata');
            continue
        end
        
        %passed all tests
        Metadata.mesh.meshing_succeeded = true;
        succeeded(c) = true;
        reason(c) = 0;
        
        if succeeded(c)
            if input.performance.verbose
            disp([paths(path_ind).mesh_file,'    seems OK']);
            end
            if plot_goods 
                figure(35)
                cla
                [s,e] = plot_mesh(Mesh,nrefines);
                %  set(s,'facealpha',1);
                title({paths(path_ind).mesh_file},'interpreter','none');
                drawnow
                pause
            end
        end
        
        
        save(paths(path_ind).metafile,'Metadata');
        
    end
    
end