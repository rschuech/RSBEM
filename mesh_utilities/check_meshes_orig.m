%%

nrefines = 1;
%false alarms:  8.6 .15     9.4 .15    10  0.85
%close to sneaking past:  6.1  0.8    6.8 0.8    5.5 0.85    9.8 0.95
volume_tol = 0.005;  % <= 0.5% error allowable
%area_outlier_tol = 10;
max_angle = pi/2; %max angle between adjacent element normals, located at centroids of the elements

if strcmp(geom(1).shape,'curved_rod') || strcmp(geom(1).shape,'tail')  %need to check each curved rod mesh for self intersection
    % disp('Waiting for Salome to finish');
    %pause;  %wait until Salome finishes everything
    
  succeeded = false(1,length(geom));
    reason = NaN(1,length(geom));
    
    
    %bad:  1516 till....1529
    
    % 0.82143 success
    
    for c = 1:length(geom)
        
%                 if succeeded(c)
%                     continue
%                 end
       
        c/length(geom)
        
        path_ind = c;  %must be c for shat to work since path variable contains full path including name
        
        
      % namebase = [outputprefix,'_','AR1','_',num2str(geom(c).AR1),'_','AR2','_',num2str(geom(c).AR2)];
        
        
        %meshfile = [paths(c).outfolder,namebase,suffix,'.dat']; %mesh dat file that comes out of Salome
        meshname = paths(path_ind).mesh_file;
        [Mesh, Metadata] = load_mesh(paths(path_ind).mesh_file);  %all paths should be the same for a given shape
        if isempty(Metadata)
            disp('empty Metadata')
            %pause
        else
            if isfield(Metadata.mesh, 'meshing_succeeded') && Metadata.mesh.meshing_succeeded
                succeeded(c) = true;
                continue
            end
        end
        
        
        
        
        if isempty(Mesh)  % || Mesh.n_elem > 1600
            disp([paths(path_ind).mesh_file,'      never created']);
           % pause
            Metadata.mesh.meshing_succeeded = false;
            save(paths(c).metafile,'Metadata');
            succeeded(c) = false;
            reason(c) = -1;
            % pause
            
            continue
            
        end
        
           
        %Metadata should always exist, maybe not Mesh
        [Mesh] = global_inds(Mesh);
        [Mesh] = renumber_Mesh(Mesh);
        clear temp_input
        temp_input.performance.nthreads = 8;
        temp_input.accuracy.integration_tol.area.abstol = 1E-9 ;
        temp_input.accuracy.integration_tol.area.reltol = 0.1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        temp_input.accuracy.integration_tol.area.maxevals = Inf;
        temp_input.accuracy.integration_tol.centroid.abstol = 1E-9 ;
        temp_input.accuracy.integration_tol.centroid.reltol = 0.1; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        temp_input.accuracy.integration_tol.centroid.maxevals = Inf;
        temp_input.accuracy.integration_tol.volume.abstol = 1E-9 ;
        temp_input.accuracy.integration_tol.volume.reltol = 0.1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        temp_input.accuracy.integration_tol.volume.maxevals = Inf;
        [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
        
        
        
     
        
        %         if isfield(Metadata.mesh,'topology')
        %             disp([paths(c).mesh_file,'     already interrogated; skipping']);
        %             if strcmp(Metadata.mesh.topology,'sphere')
        %                 succeeded(c) = true;
        %             else
        %                 succeeded(c) = false;
        %             end
        %             continue
        %         end
        
        
        
        edges = get_edges(Mesh);  %used for both topology check as well as adjacent normals check
        
        metric = size(unique(Mesh.elems(:,1:3)),1) - size(edges.verts,1) + Mesh.n_elem;  % # corner verts - # edges + # elements     = 2 for sphere, 0 for torus
        
        if metric ~= 2 %don't have a topological sphere, something is wrong, probably self-intersection
            if metric == 0
                disp([paths(path_ind).mesh_file,'    is a topological torus']);
            else
                disp([paths(path_ind).mesh_file,'    is topologically wierd']);
            end
            
            figure(34)
            cla
            [s,e] = plot_mesh(Mesh,nrefines);
            set(s,'facealpha',1);
            title(paths(path_ind).mesh_file,'interpreter','none');
            drawnow
            %   disp('Mark mesh as bad and continue?');
           % pause
            
            
            
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
            %  Metadata.geom.self_intersected = true;
            
            %keep Mesh file around in case you want to look at it later
            %             %delete bad mesh file
            %             delete(paths(c).mesh_file);
        else %mesh appears OK
            %Metadata.geom.self_intersected = false;
            
            Metadata.mesh.topology = 'sphere';
            %wait, mesh still may be fucked up in other ways
            %             Metadata.mesh.meshing_succeeded = true;
            %             succeeded(c) = true;
            %             reason(c) = 0;
        end
        
        % not out of the woods yet; may still have the wrong volume
        if ~strcmp(geom(1).shape,'tail') && (  Mesh.Volume < (1 - volume_tol) || Mesh.Volume > (1 + volume_tol)  ) %if volume is either too small or too large
            disp([paths(path_ind).mesh_file,'    has volume problem, error = ',num2str(100*(1-Mesh.Volume)),' %']);
            figure(34)
            cla
            [s,e] = plot_mesh(Mesh,nrefines);
            set(s,'facealpha',1);
            title({paths(path_ind).mesh_file,['has volume problem, error = ',num2str(100*(1-Mesh.Volume)),' %']},'interpreter','none');
            drawnow
           % pause
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
              disp([paths(path_ind).mesh_file,'    has normals angle problem']);
                figure(34)
            cla
            [s,e] = plot_mesh(Mesh,nrefines);
            set(s,'facealpha',0.6);
            title({paths(path_ind).mesh_file,['has normals angle problem']},'interpreter','none');
            drawnow
%             pause
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
            disp([paths(path_ind).mesh_file,'    seems OK']);
          
%                              figure(34)
%                         cla
%                         [s,e] = plot_mesh(Mesh,nrefines);
%                        %  set(s,'facealpha',1);
%                         title({paths(path_ind).mesh_file},'interpreter','none');
%                         drawnow
%                       pause
        end
        
        
        save(paths(path_ind).metafile,'Metadata');
        
    end
    
end