

      [field_vel] = field_velocity_mexed(Mesh_clearance.verts, f, Mesh, field_input);  % Mesh is the full mesh correspondng to solution f
        % points = N x 3
        % field_vel = N x 3




function [ ] = 



field_input.accuracy.integration_tol.field_vel.abstol = 1E-7;
field_input.accuracy.integration_tol.field_vel.abstol = 1E-6;
%field_input.accuracy.integration_tol.field_vel.abstol = 1E-5;

field_input.accuracy.integration_tol.field_vel.reltol = 0;
field_input.accuracy.integration_tol.field_vel.maxevals = Inf;
field_input.accuracy.eps2 = input.accuracy.eps2;
field_input.constants.multfactor = input.constants.multfactor;
field_input.constants.mu = input.constants.mu;

folder = 'C:\Users\rudi\Desktop\RD\dino_dumps_hairs_3.5_1.5_phases\';

dump_folder = 'C:\Users\rudi\Desktop\RD\base_case_2_1.5\';
mesh_folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\regular\final\';

outfolder = 'C:\Users\rudi\Desktop\RD\dino velocity fields\hairs 2 1.5\';


dumps = dir([dump_folder,'*dump.mat']);
dumps = {dumps.name};

for d = 1 %1:length(dumps)
    [d length(dumps)];
    dump = dumps{d}
    
    
    load([dump_folder,dump],'Mesh_files','Solutions','input');
    if exist('mesh_folder','var') && ~isempty(mesh_folder)
        input.paths.datfolder = mesh_folder;
    end
    
    clear vel_fields

    
    % has Solutions struct
    Solutions_fields = setdiff( fieldnames(Mesh_files) , {'phase','time','metadata'});  %order of submeshes as far as Solutions (e.g. traction f) is concerned
    % hmm, it actually appears that the order of Solutions_fields is the same
    % as the order of Solutions.rand_inds, so probably don't need to do
    % anything to match the orders
    
    %%
    Xlist = [];  Ylist = [];  Zlist = [];
    
  
    
    points = [Xlist(:) Ylist(:) Zlist(:)];
    points = unique(points,'rows');
    
    Meshes = [];
    
    % average over a phase cycle
    fh = waitbar(0,'phase cycle velocity field avg');
    for phase_ind = 1:length(Solutions.phase)
        waitbar(phase_ind / length(Solutions.phase) );
        disp(['On phase ind ',num2str(phase_ind),' of ',num2str(length(Solutions.phase))]);
    
        phase = Solutions.phase(phase_ind);
        
        [Mesh, Metadata, matrix_props, Metadata_rand_inds] = load_dino_mesh(phase, Mesh_files, input, Solutions.rand_inds{phase_ind});
        
        if do_rotation
            [Mesh,rotmat] = rotateMesh(Mesh, rotation_angle);
            % Metadata apparently never used here so don't bother rotating
            % it
          
            Solutions.f{phase_ind} = reshape( rotmat * reshape( Solutions.f{phase_ind} , 3, [])  , [] , 1);  % reshape tracton vector to a matrix, rotate it, reshape back to a vector
            
            Solutions.U(phase_ind,:) = (rotmat * Solutions.U(phase_ind,:)')';
            Solutions.Omega(phase_ind,:) = (rotmat * Solutions.Omega(phase_ind,:)')';
        end
        
        Meshes{phase_ind} = Mesh;
    
        
        tic;
        
        
        % Lasse tested 1 - 30 um dia prey particles in mixotrophic pape
        % assume prey captured when it comes within 1 prey radius of dino
        % surface
        
        % points = N x 3
        % field_vel = N x 3
        [field_vel] = field_velocity_mexed(points, Solutions.f{phase_ind}, Mesh, field_input);
        
        
        
        
        disp(['field_velocity_mexed took ',num2str(toc)]);
        
        U = NaN(size(X));  V = U;  W = U;
        parfor i = 1:numel(X)
            
            point = [X(i) Y(i) Z(i)];  %point in full 3D grid
            ind = find(  points(:,1) == point(1) & points(:,2) == point(2) & points(:,3) == point(3) );
            if ~isempty(ind)
                if numel(ind) > 1
                    stopafra
                end
                
                U(i) = field_vel(ind,1);  V(i) = field_vel(ind,2);  W(i) = field_vel(ind,3);
                %    pause
            else
                U(i) = NaN;  V(i) = NaN;  W(i) = NaN;  %no idea why this is needed - parfor is forgetting the initialization of U,V,W and defaulting it to zeros!
            end
        end
        
        vel_fields(phase_ind).phase = Solutions.phase(phase_ind);
        vel_fields(phase_ind).X = X; vel_fields(phase_ind).Y = Y;  vel_fields(phase_ind).Z = Z;
        vel_fields(phase_ind).U = U; vel_fields(phase_ind).V = V;  vel_fields(phase_ind).W = W;
        vel_fields(phase_ind).field_vel = field_vel;
        vel_fields(phase_ind).points = points;
        
    end
    close(fh)
    
    
    %%

    names = {'U','V','W'};
    vel_field_avg_U = NaN(size(U));  vel_field_avg_V = NaN(size(V));  vel_field_avg_W = NaN(size(W));

    parfor i = 1:numel(vel_fields(1).X)
        temp = []; temp2 = []; fun = [];  fun2 = [];
        for n = 1:length(names)
            j = [];
            y = []; yb = [];
            for j = 1:length(Solutions.phase)
                y(j) = vel_fields(j).(names{n})(i);
                yb(j) = vel_fields(j).bodyframe.(names{n})(i);
            end
            switch input.kinematics_interp_method
                case 'trig'
                    temp = trig_interp_fit(Solutions.phase,y);
                    temp2 = trig_interp_fit(Solutions.phase,yb);

                case 'spline'
                    
                    filter = ~ismember(roundn(Solutions.phase,input.interp_phases_tol), roundn(input.discard_interp_phases,input.interp_phases_tol));
                   if ~all(isnan(y(filter)))
                    temp = csape([Solutions.phase(filter) ; phase_bounds(2)]' , [y(filter)  y(1)],'periodic'); % assume that first entry is phase = 0, which is never discarded
                    temp2 = csape([Solutions.phase(filter) ; phase_bounds(2)]' , [yb(filter)  yb(1)],'periodic');
                   else
                       temp = [];  temp2 = [];
                   end

                    
            end
            
            
            if ( strcmp( input.kinematics_interp_method , 'trig') && isnan(temp.a0) ) || (  strcmp( input.kinematics_interp_method , 'spline') && isempty(temp)  )
                vel_field_avg_U(i) = NaN;  vel_field_avg_V(i) = NaN;  vel_field_avg_W(i) = NaN;
                vel_field_avg_U_bodyframe(i) = NaN;  vel_field_avg_V_bodyframe(i) = NaN;  vel_field_avg_W_bodyframe(i) = NaN;
            else
                switch input.kinematics_interp_method
                    case 'trig'
                        fun = @(x) trig_interp_eval(temp, x);
                        fun2 = @(x) trig_interp_eval(temp2, x);
                    case 'spline'
                        fun = @(x) fnval(temp,x);
                        fun2 = @(x) fnval(temp2,x);
                end
                
                switch names{n}
                    case 'U'
                        vel_field_avg_U(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                        vel_field_avg_U_bodyframe(i) = 1/diff([0 2*pi]) * integral(fun2, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                    case 'V'
                        vel_field_avg_V(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                        vel_field_avg_V_bodyframe(i) = 1/diff([0 2*pi]) * integral(fun2, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                    case 'W'  
                        vel_field_avg_W(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                        vel_field_avg_W_bodyframe(i) = 1/diff([0 2*pi]) * integral(fun2, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                end
            end
        end
        
    end
    
    save([outfolder,dumps{d}(1:end-8),'vel_fields_z_slice.mat']);
    return
    %%
    min_distance = 2;  do_waitbar = true;
    [OK_points] = filter_inside_or_near_boundaries(points,min_distance,Solutions,Meshes, do_waitbar);
    %%
end

return

%%


        
 