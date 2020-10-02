function [RHS] = assemble_RHS(input,Mesh, matrix_props, BCs)

RHS = zeros(matrix_props.n_rows,1);

rhstic = tic;
switch input.problemtype
    case 'forced'
        vertlist = [];  normals = [];
        for i = 1:length(Mesh)
            inds = find( Mesh(i).indices.glob.vert >= Mesh(i).indices.glob.unq_bounds.vert(1) &   Mesh(i).indices.glob.vert <= Mesh(i).indices.glob.unq_bounds.vert(2)); %local inds of "global" verts belonging to this mesh
            vertlist = [vertlist; Mesh(i).verts(inds,:)];   %methinks global indices are always in increasing order, so this always works (?)
            normals = [normals;  Mesh(i).normals(inds,:)];
            
        end
        
       
        for i_mesh_vert = 1:length(Mesh)
            col_inds = Mesh(i_mesh_vert).indices.glob.unq_bounds.vert(1) : Mesh(i_mesh_vert).indices.glob.unq_bounds.vert(2);
            switch input.BC_type(i_mesh_vert)
                case 1 % no slip
                    temp = velocities(col_inds,:);
                    temp = temp';  %switch to rows being x, y, z and cols being pts
                    RHS(3*(col_inds(1)-1)+1 : 3*(col_inds(end)-1)+3) = temp(:);  %works since each column is placed under the one before, and each column is [u v w]' for a pt
                    
                    
                case 2 % free slip
                    
                  %  velocities = repmat([0* 50/sqrt(3)],size(velocities,1),3);
                    RHS(3*(col_inds-1)+1) = velocities(col_inds,1) .* normals(col_inds,1) + velocities(col_inds,2) .* normals(col_inds,2) + velocities(col_inds,3) .* normals(col_inds,3);
                    % no need to do anything for f.s1 = 0, f.s2 = 0 since RHS stays
                    % zero for those
                    
            end
            
        end
        
        
        
    case 'freeswim'
        
        %RHS(1:(Mesh(1).n_vert)*3,1) = 0; %body traction RHS is always 0 since U and Omega are always unknown
        %above line unnecessary due to initialization with zeros
        switch input.bugtype
            case 'bacteria'
                switch input.tail.motorBC
                    case 'freq'  %freq is known, goes into RHS for tail traction
                        temp = zeros(Mesh(2).indices.glob.unq_bounds.vert(2) - Mesh(2).indices.glob.unq_bounds.vert(1) + 1,3);
                        r = Mesh(2).verts - repmat(Mesh(2).refpoints(:,1)',Mesh(2).n_vert,1);  %use tail refpoint since it's always on motor axis
                        temp(:,1) = BCs.freeswim.motor_freq*(Mesh(1).orientation(2)*r(:,3) - Mesh(1).orientation(3)*r(:,2));  %x equations
                        temp(:,2) = BCs.freeswim.motor_freq*(Mesh(1).orientation(3)*r(:,1) - Mesh(1).orientation(1)*r(:,3));  %y equations
                        temp(:,3) = BCs.freeswim.motor_freq*(Mesh(1).orientation(1)*r(:,2) - Mesh(1).orientation(2)*r(:,1));  %z equations
                        
                        % don't seem to currently allow body mesh to be
                        % free slip, would need to do that here....?
                        switch input.BC_type(2)
                            case 1 % no slip
                                RHS((Mesh(2).indices.glob.unq_bounds.vert(1)-1)*3+1  :3:  Mesh(2).indices.glob.unq_bounds.vert(2)*3,1) = temp(:,1);
                                RHS((Mesh(2).indices.glob.unq_bounds.vert(1)-1)*3+2  :3:  Mesh(2).indices.glob.unq_bounds.vert(2)*3,1) = temp(:,2);
                                RHS((Mesh(2).indices.glob.unq_bounds.vert(1)-1)*3+3  :3:  Mesh(2).indices.glob.unq_bounds.vert(2)*3,1) = temp(:,3);
                                
                                %for dino, just need to change above to take known surface
                                %velocities in orig body frame, and rotate them just like we
                                %rotated the orientation vector for the body in yprime.m or
                                %whatever
                                
                            case 2
                                inds = find( Mesh(2).indices.glob.vert >= Mesh(2).indices.glob.unq_bounds.vert(1) &   Mesh(2).indices.glob.vert <= Mesh(2).indices.glob.unq_bounds.vert(2)); %local inds of "global" verts belonging to this mesh
                                normals = Mesh(2).normals(inds,:);
                                RHS((Mesh(2).indices.glob.unq_bounds.vert(1)-1)*3+1  :3:  Mesh(2).indices.glob.unq_bounds.vert(2)*3,1) = temp(:,1).*normals(:,1) + temp(:,2).*normals(:,2) + temp(:,3).*normals(:,3);
                        end
                        
                    case 'torque'  %omega is unknown, RHS is almost all 0 since freq is unknown but have additional equation for last row
                        %RHS(Mesh(2).global_indices.vert.start : Mesh(2).global_indices.vert.end*3,1) = 0;
                        %above unnecessary due to zero initialization of RHS
                        RHS(Mesh(2).indices.glob.unq_bounds.vert(2)*3+7,1) = BCs.freeswim.motor_torque;
                end
            case {'dino', 'sheet'}
                % velocities are stored in BCs, which is rotated along with
                % Mesh inside yprime.m, so here it is already in the same
                % orientation as Mesh and can be added right in to RHS
                
                
                % tot_verts = Mesh(end).indices.glob.unq_bounds.vert(2 );
                velocities = [] ;  %x, y, z velocity at each unq vert
                normals = [];
                for i = 1:length(Mesh)
                    inds = find( Mesh(i).indices.glob.vert >= Mesh(i).indices.glob.unq_bounds.vert(1) &   Mesh(i).indices.glob.vert <= Mesh(i).indices.glob.unq_bounds.vert(2)); %local inds of "global" verts belonging to this mesh
                    velocities = [velocities; BCs.freeswim.(Mesh(i).name)(inds,:)];
                    normals = [normals;  Mesh(i).normals(inds,:)];
                end
                
                for i_mesh_vert = 1:length(Mesh)
                    col_inds = Mesh(i_mesh_vert).indices.glob.unq_bounds.vert(1) : Mesh(i_mesh_vert).indices.glob.unq_bounds.vert(2);
                    switch input.BC_type(i_mesh_vert)
                        case 1 % no slip
                            temp = velocities(col_inds,:);
                            temp = temp';  %switch to rows being x, y, z and cols being pts
                            RHS(3*((col_inds(1))-1)+1 : 3*((col_inds(end))-1)+3) = temp(:);  %works since each column is placed under the one before, and each column is [u v w]' for a pt
                            
                            
                            
                        case 2 % free slip
                            RHS(3*((col_inds)-1)+1) = velocities(col_inds,1) .* normals(col_inds,1) + velocities(col_inds,2) .* normals(col_inds,2) + velocities(col_inds,3) .* normals(col_inds,3);
                            % no need to do anything for f.s1 = 0, f.s2 = 0 since RHS stays
                            % zero for those
                            
                            
                            %last 6 rows are always zero (as initialized) since they correspond to
                            %sum(forces) = 0 and sum(torques) = 0 for any freeswim
                            %problem
                    end
                    
                end
                
        end
        
end

% disp(['RHS assembly took ',num2str(toc(rhstic))]);  %takes almost no time