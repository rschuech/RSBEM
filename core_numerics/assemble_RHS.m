function [RHS] = assemble_RHS(input,Mesh, matrix_props, BCs)


rhstic = tic;
switch input.problemtype
    case 'forced' %for now, entire object is treated as rigidly connected submeshes
        vertlist = [];
        for i = 1:length(Mesh)
            inds = find( Mesh(i).indices.glob.vert >= Mesh(i).indices.glob.unq_bounds.vert(1) &   Mesh(i).indices.glob.vert <= Mesh(i).indices.glob.unq_bounds.vert(2)); %local inds of "global" verts belonging to this mesh
            vertlist = [vertlist; Mesh(i).verts(inds,:)];   %methinks global indices are always in increasing order, so this always works (?)
        end
        
        %forced rotations will occur around refpoint
        if strcmp(input.bugtype,'bacteria') && length(Mesh) == 2
            refpoint = Mesh(2).refpoints(:,1); %make sure refpoint is somewhere along the axis of motor rotation
        else
            refpoint = Mesh(1).refpoints(:,1); % doesn't really matter where it is, so use body or whatever first submesh is
        end
        
        [u] = compute_BCs_mexed(BCs.forced,vertlist,refpoint,input.performance.nthreads);
        RHS = u;  %for forced movement, RHS of matrix equation equals u
    case 'freeswim'
        RHS = zeros(matrix_props.n_rows,1); %initialize
        
        %RHS(1:(Mesh(1).n_vert)*3,1) = 0; %body traction RHS is always 0 since U and Omega are always unknown
        %above line unnecessary due to initialization with zeros
        switch input.bugtype
            case 'bacteria'
                switch input.tail.motorBC
                    case 'freq'  %freq is known, goes into RHS for tail traction
                        r = Mesh(2).verts - repmat(Mesh(2).refpoints(:,1)',Mesh(2).n_vert,1);  %use tail refpoint since it's always on motor axis
                        RHS((Mesh(2).indices.glob.unq_bounds.vert(1)-1)*3+1  :3:  Mesh(2).indices.glob.unq_bounds.vert(2)*3,1) = BCs.freeswim.motor_freq*(Mesh(1).orientation(2)*r(:,3) - Mesh(1).orientation(3)*r(:,2));  %x equations
                        RHS((Mesh(2).indices.glob.unq_bounds.vert(1)-1)*3+2  :3:  Mesh(2).indices.glob.unq_bounds.vert(2)*3,1) = BCs.freeswim.motor_freq*(Mesh(1).orientation(3)*r(:,1) - Mesh(1).orientation(1)*r(:,3));  %y equations
                        RHS((Mesh(2).indices.glob.unq_bounds.vert(1)-1)*3+3  :3:  Mesh(2).indices.glob.unq_bounds.vert(2)*3,1) = BCs.freeswim.motor_freq*(Mesh(1).orientation(1)*r(:,2) - Mesh(1).orientation(2)*r(:,1));  %z equations
                        
                        %for dino, just need to change above to take known surface
                        %velocities in orig body frame, and rotate them just like we
                        %rotated the orientation vector for the body in yprime.m or
                        %whatever
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
                temp = [] ;  %x, y, z velocity at each unq vert
                for i = 1:length(Mesh)
                    inds = find( Mesh(i).indices.glob.vert >= Mesh(i).indices.glob.unq_bounds.vert(1) &   Mesh(i).indices.glob.vert <= Mesh(i).indices.glob.unq_bounds.vert(2)); %local inds of "global" verts belonging to this mesh
                    temp = [temp; BCs.freeswim.(Mesh(i).name)(inds,:)];
                end
                
                RHS(1:3:matrix_props.n_col*3) = temp(:,1);  %slot in x-velocities
                RHS(2:3:matrix_props.n_col*3) = temp(:,2);  %slot in y-velocities
                RHS(3:3:matrix_props.n_col*3) = temp(:,3);  %slot in z-velocities
                %last 6 rows are always zero (as initialized) since they correspond to
                %sum(forces) = 0 and sum(torques) = 0 for any freeswim
                %problem
    
        end
        
end

% disp(['RHS assembly took ',num2str(toc(rhstic))]);  %takes almost no time