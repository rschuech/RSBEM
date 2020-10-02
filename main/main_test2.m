% Inputs = Inputs(1:117);

%clearvars -except Inputs
inputs_loop_tic = tic;
%if run from Condor, current directory will be temporary folder where
%Condor has placed wrapper.m
%cd E:/Hull/code/;
%cd C:\RD\bem\code\;
%addpath(genpath('./')); %add all subfolders to path so that subfunctions will be found
solutions_temp = [];
D = [];

%bad_sweep_i = false(1,length(Inputs));  %keep track of which runs were bad due to not having mesh files
for sweep_i = 1:length(Inputs)
    %     if sweep_i == 139
    %         continue
    %     end
    
    
    input = Inputs(sweep_i);
    if input.performance.verbose
        disp(['On sweep iter ',num2str(sweep_i),' of ',num2str(length(Inputs)),'']);
    end
    total_time_tic = tic;
    
    clear Mesh Metadata index_mapping  timings BCs
    
    %% all this is too ensure we can reuse part of A matrix between resistance and mobility, in conjunction with randomize_verts
    if sweep_i > 1  %can't reuse anything if this is the first sweep iteration, duh
        [diff_fields] = struct_diff(Inputs(sweep_i-1),Inputs(sweep_i));  %if we can reuse, prev inputs should be identical except for a few things:
        input_diffs = setdiff(diff_fields,{'full','problemtype'});
    else
        input_diffs = NaN;  %anything that's not empty works
    end
    
    if isempty(input_diffs)  %we just did resistance and are now doing mobility or vice versa
        if is_intersected  %can use same results of intersections check as from last run, since it was same geometry
            disp('Tail intersects Body (known from previous run); aborting.');
            continue
        end
        recycle_A = true;
    else
        recycle_A = false;
    end
    
    switch input.bugtype
        case {'bacteria'}
            if ~exist('persistant_data','var')
                persistant_data(1).rand_inds = [];  %use placeholder to start; this will hold rand_inds for each submesh and maybe other things to save between sweep iters later
                persistant_data(2).rand_inds = [];
            end
        case {'dino', 'sheet'}
            persistant_data(1).rand_inds = [];  %use placeholder to start; this will hold rand_inds for each submesh and maybe other things to save between sweep iters later
            persistant_data(2).rand_inds = [];
            persistant_data(3).rand_inds = [];
            persistant_data(4).rand_inds = [];
            persistant_data(5).rand_inds = [];
            persistant_data(6).rand_inds = [];
    end
    %%
    
    %% load body mesh
    
    
    switch input.bugtype
        case {'bacteria'}
            if strcmp(input.problemtype,'mobility') || ( strcmp(input.problemtype,'resistance') && input.include_tail )
                submeshes = ["Body" , "Tail"];  %which submeshes to load  (always list body first by convention)
                %                 rigid = [true true]; % this refers to whether the *mesh* deforms over time (note that geometry can be rigid but mesh may still
                % deform/change across time, e.g. dino body, since it currently gets remeshed along with transverse sheet at every phase)
            else
                submeshes = "Body";
                %                 rigid = true;
            end
        case 'dino'
            % temp = {'Body','Transverse','Tail','Wingtip'};
            %             submeshes = {'body','transverse','tail'};
            %  submeshes = temp(logical(input.potatohead));  %only load submeshes that will be simulated - not sure if this will break shat for old cases.....
            
            % changed to automatically figure out which submeshes were used
            % during compute BCs internally, and load all of them
            [Mesh_files] = get_dino_mesh_files(input.paths.datfolder);  %Mesh_files is a list of the mesh and Metadata filenames for all the phase angles
        case "sheet"
            submeshes = "sheet";
            [Mesh_files] = get_sheet_mesh_files(input.paths.datfolder);
            %             rigid = false;
    end
    
    
    %for monoflagellated bacteria model, Mesh(1) is the body, Mesh(2) is the tail
    if strcmp(input.bugtype,'bacteria')   || (  strcmp(input.bugtype,'dino') && strcmp(input.problemtype,'resistance')  )  %dino mobility case has meshes loaded again later anyway in interp_y
        bad_sweep_i(sweep_i) = false;
        for si = 1:length(submeshes)
            meshname = input.paths.datfolder + input.paths.namebase.(submeshes(si)) + '.dat';
            if input.performance.verbose
                disp(['Loading ' + submeshes(si) + ' mesh     ' + meshname]);
            end
            
            
            try
                %                 [temp_mesh, Metadata(si), persistant_data(si).rand_inds] = load_mesh_wrapper(meshname, input, persistant_data(si).rand_inds, input_diffs);  %old rand_inds goes in, new rand_inds comes out (which might be the same as it was)
                [temp_mesh, temp_metadata,  index_mapping.local_node2global_node{si}] = load_mesh(meshname, input);
                if ~isempty(temp_metadata)
                    Metadata(si) = temp_metadata;
                end
            catch
                [mesh_succeed] = gen_mesh_wrapper(input);
                if ~mesh_succeed
                    stopafra
                end
                [temp_mesh, Metadata(si), index_mapping.local_node2global_node{si}] = load_mesh(meshname, input);
                %                 [temp_mesh, Metadata(si), persistant_data(si).rand_inds] = load_mesh_wrapper(meshname, input, persistant_data(si).rand_inds, input_diffs);  %old rand_inds goes in, new rand_inds comes out (which might be the same as it was)
            end
            
            %             if ( isfield(Metadata(si).mesh,'meshing_succeeded') && ~Metadata(si).mesh.meshing_succeeded ) || isempty(temp_mesh)
            %                 bad_sweep_i(sweep_i) = true;
            %                 break %don't bother trying to load other submeshes
            %             end
            
            temp_mesh.name = submeshes(si);
            temp_mesh.parent_topology = input.parent_topology.( temp_mesh.name );
            temp_mesh.is_mesh_rigid = input.is_mesh_rigid.(temp_mesh.name);
            temp_mesh.orientation = [1 0 0; 0 1 0; 0 0 1;]';  %initial orientation vectors are along x, y, z axes
            
            Mesh(si) = temp_mesh;  clear temp_mesh
            
        end
        
        if bad_sweep_i(sweep_i)
            continue
        end
        
    end
    
    
    Mesh0 = Mesh;
    % eliminate any possible overlap in global indices across coincident submesh groups, and also make sure global indices are consecutive
    [Mesh, index_mapping.local_node2global_node] = shift_global_indices(Mesh, index_mapping.local_node2global_node, input.coincident_submeshes);
    
    if strcmp(input.bugtype,'bacteria')  || (  strcmp(input.bugtype,'dino') && strcmp(input.problemtype,'resistance')  )
        %         [Mesh] = global_inds(Mesh); % adds global values to Mesh.indices
        %         [Mesh] = renumber_Mesh(Mesh);
        
        % create smaller input struct with just necessary parameters to reduce mex
        % pain
        %         temp_input.performance.nthreads = input.performance.nthreads;
        %         temp_input.accuracy.integration_tol.area = input.accuracy.integration_tol.area;
        %         temp_input.accuracy.integration_tol.centroid = input.accuracy.integration_tol.centroid;
        %         temp_input.accuracy.integration_tol.volume = input.accuracy.integration_tol.volume;
        %         [Mesh] = store_mesh_cons+tants_wrapper(Mesh, temp_input);
        
        [Mesh] = compute_element_based_mesh_parameters(Mesh, input);
        [node_parameters, index_mapping, Mesh] = compute_node_based_parameters(Mesh,index_mapping, input.parent_topology, input.coincident_submeshes, input.BC_type);
    end
    
    
    
    if strcmp(input.bugtype,'bacteria') && strcmp(input.Body.shape,'curved_rod')
        
        if input.Body.AR(1) == 1 %sphere degenerate case
            Mesh(1).refpoints = [0 0 0]';   %center
        elseif input.Body.AR(2) == 0  %not a sphere, but a straight rod degenerate case
            Mesh(1).refpoints = [0 0 0]';   %center
        elseif isfield(Metadata,"geom")  %actually a general curved rod
            Mesh(1).refpoints =  [2*Metadata(1).geom.radius_curv*sin(Metadata(1).geom.nturns / 2 *pi) * cos(Metadata(1).geom.nturns / 2 *pi); ...
                2*Metadata(1).geom.radius_curv*(sin(Metadata(1).geom.nturns / 2 *pi))^2;
                0 ];   %center of centerline
        else
            Mesh(1).refpoints = [0 0 0]';   %probably center
        end
    else
        Mesh(1).refpoints = [0 0 0]';   %probably center
    end
    
    if length(Mesh) == 2 && input.include_tail
        Mesh(2).refpoints = [0 0 0]';  % will be used for torque calculations - OK as long as somewhere along motor axis
    end
    
    
    
    if  strcmp(input.bugtype,'bacteria') && input.include_tail && input.Tail.rotate_tail %to check for possible effect of tail phase angle on calculated diffusivities   %strcmp(input.problemtype,'resistance') &&
        % do this even for mobility problem so that if the tail was rotated for
        % resistance, it will get rotated for mobility problem and thus part of A
        % matrix can be reused.  Not that we're likely to ever do
        % mobility problem runs in conjunction with this option anyway....
        [Mesh(2)] = rotateMesh(Mesh(2), [   0  0  input.tail.rotate_tail_angle]' );  %rotate tail around x axis
    end
    
    if strcmp(input.bugtype,'bacteria')
        switch input.Tail.orientation
            case 'pole2pole'
                [~,sphere_centers, centerline_center] = calc_pole_coords(Metadata);
                
                if Metadata(1).geom.AR1 == 1 %sphere degenerate case, origin is body center
                    rear_to_front_dir = [1 0 0]'; % need to hardcode since head and tail sphere centers are same point
                else
                    rear_to_front_dir = sphere_centers.head - sphere_centers.tail;  rear_to_front_dir = rear_to_front_dir / sqrt(sum(rear_to_front_dir.^2));
                end
                pole2pole_line = [sphere_centers.tail' rear_to_front_dir'];
                
                if Metadata(1).geom.AR1 == 1 %sphere degenerate case, origin is body center
                    rear_pole_circle = [0 0 0 Metadata(1).geom.radius];
                elseif Metadata(1).geom.AR2 == 0  %not a sphere, but a straight rod degenerate case, origin is body center
                    
                    rear_pole_circle = [-Metadata(1).geom.height/2  0 0 Metadata(1).geom.radius];
                else  %actually a general curved rod, origin is center of tail-end sphere
                    rear_pole_circle = [0 0 0 Metadata(1).geom.radius];
                    
                end
                
                point = intersectLineSphere(pole2pole_line, rear_pole_circle);
                point( point(:,1) > rear_pole_circle(1) | point(:,2) > rear_pole_circle(2), :) = NaN;  % any intersection points not within lower left circle quadrant are not of interest
                point( all(isnan(point),2) , :) = []; % remove any non-interesting intersection points
                
                if isempty(point) % must instead intersect along main body segment
                    centerline_circle = [centerline_center'  (Metadata(1).geom.radius_curv + Metadata(1).geom.radius)];
                    point = intersectLineSphere(pole2pole_line, centerline_circle);
                    point( point(:,1) < centerline_circle(1) | point(:,2) > rear_pole_circle(2), :) = NaN;
                    point( all(isnan(point),2) , :) = [];
                end
                
                
                if isempty(point)
                    error('Couldn''t find intersection point');
                end
                
                if size(point,1) > 1
                    error('More than one feasible intersection point found.');
                end
                
                angle = acos( dot(Mesh(2).orientation(:,1) , rear_to_front_dir) );
                
                
                
                tail_proximal_pole = Mesh(2).refpoints + [Metadata(2).geom.pipeRadius 0 0]';
                Mesh(2) = shiftMesh(Mesh(2),[ - tail_proximal_pole]); % shift so that proximal pole point is at origin in order to rotate easily
                Mesh(2) = rotateMesh(Mesh(2), angle, [0 0 1]');
                
                Mesh(2) = shiftMesh(Mesh(2),[point' - 0]); % proximal pole point is the origin right now
                
        end
    end
    
    %     angle.vertical = 0.5  * -pi/2;  %"declination angle" should be between 0 and -pi/2 since positive would be same due to symmetry, and more than 90 deg is slightly ridiculous
    %     angle.horizontal = 0.5 * -pi/2;  %"azimuthal angle" should be between -pi/2 and pi/2, more than this would be repeating due to symmetry
    %
    %     [Mesh(2)] = rotateMesh(Mesh(2),[angle.horizontal, angle.vertical, 0],[2 1 3] );
    if strcmp(input.bugtype,'bacteria') && length(Mesh) == 2
        tail_mesh = Mesh(2);
        tail_shift = 0;
        tail_shift_increment = input.accuracy.check_intersections_tolerance;
    end
    break_switch = false;
    while ~break_switch
        
        %%
        if strcmp(input.bugtype,'bacteria') && length(Mesh) == 2
            %perform translations to get tail next to body
            %start tail one tail radius away from body
            switch input.Body.shape
                case 'ellipsoid' %origin for body located at body center
                 
                    Mesh(2) = shiftMesh(Mesh(2),[-2*input.Tail.radius - Metadata(1).geom.a 0 0]);
                
                case 'capsule' %origin for body located at body center
                    Mesh(2) = shiftMesh(Mesh(2),[-2*input.tail.radius - (Metadata(1).geom.caps_c + Metadata(1).geom.height/2) 0 0]);
                case 'curved_rod'  %origin for body located at center of spherical cap
                    switch input.Tail.orientation
                        case 'centerline'
                            if input.Body.AR(1) == 1 %sphere degenerate case
                                Mesh(2) = shiftMesh(Mesh(2),[-2*input.Tail.radius - Metadata(1).geom.radius 0 0]);
                            elseif input.body.AR(2) == 0  %not a sphere, but a straight rod degenerate case
                                Mesh(2) = shiftMesh(Mesh(2),[-2*input.tail.radius - (Metadata(1).geom.radius + Metadata(1).geom.height/2) 0 0]);
                            else  %actually a general curved rod
                                Mesh(2) = shiftMesh(tail_mesh,[(-2*input.tail.radius - (Metadata(1).geom.radius) - tail_shift) 0 0]);
                            end
                        case 'pole2pole'
                            dist = input.tail.radius + tail_shift;
                            Mesh(2) = shiftMesh(tail_mesh, -rear_to_front_dir * dist);
                    end
                    
            end
            
        else
            break_switch = true;
            
        end
        
        if (input.bugtype == "bacteria" && length(Mesh) > 1 )  || (  strcmp(input.bugtype,'dino') && strcmp(input.problemtype,'resistance')  )
            intersection_check_tic = tic;
            if input.Body.AR(2) <= 0.8
                %                 disp('skipping intersections check');
                is_intersected = false;
            else
                nrefines = [2 1];
                [is_intersected] = submesh_intersections(Mesh,node_parameters,index_mapping,input.accuracy.check_intersections_tolerance, input.accuracy.check_intersections_n_angles,false,input.performance.nthreads,nrefines);  % checks for self-intersection at each angle in parallel
            end
            
            timings.intersection_check = toc(intersection_check_tic);
            
            if is_intersected
                %                 bad_sweep_i(sweep_i) = true;
                %             load(input.paths.intersections_file,'intersections'); %contains one variable called intersections, a list of non-working input structs
                %             if isempty(intersections)
                %                 intersections = input;
                %             else
                %                 intersections(end+1) = input;  %should contain all the info we'll ever need on the current geometry combination that doesn't work
                %             end
                %             save(input.paths.intersections_file,'intersections');
                %                 disp('Tail intersects Body; aborting.');
                disp('Tail intersects Body; Shifting.');
                tail_shift = tail_shift + tail_shift_increment;
                %don't even bother saving anything to Results for impossible
                %configurations
                
                %                 continue %with Inputs parameter sweep loop
            else
                if input.performance.verbose
                    disp('Intersections test passed.');
                end
                break_switch = true;
            end
            
        end
        
    end  % while there are self-intersections
    
    %if computing friction coeffs, shift whole thing so that fixed ref
    %frame origin is at centroid (though this isn't really necessary)
    if (strcmp(input.problemtype,'resistance') && input.include_tail)
        Centroid_overall = sum(repmat([Mesh.Volume],3,1) .* [Mesh.Centroid],2) ./ sum([Mesh.Volume]); %Centroid of entire mesh
        Mesh = shiftMesh(Mesh, - Centroid_overall); %after this point, motor axis no longer coincides with x-axes, so don't try to rotateMesh the tail!
    end
    
    
    
    
    
    
    
    
    % make smaller version of input struct with just variables needed to assemble A
    % matrix, to reduce pain of re-mexing if input struct is modified
    clear assembly_input
    assembly_input.performance.nthreads = input.performance.nthreads;
    assembly_input.problemtype = input.problemtype;
    %     assembly_input.BC_type = input.BC_type;
    %
    %         assembly_input.BC_type = [];
    %     for field = string({Mesh.name})
    %         if isfield(input.BC_type,field)
    %             assembly_input.BC_type(end+1) = input.BC_type.(field);
    %         else
    %             assembly_input.BC_type(end+1) = input.BC_type.default;
    %         end
    %     end
    
    assembly_input.performance.eliminate_DL = logical([]);
    for field = string({Mesh.name})
        if isfield(input.performance.eliminate_DL,field)
            assembly_input.performance.eliminate_DL(end+1) = input.performance.eliminate_DL.(field);
        else
            assembly_input.performance.eliminate_DL(end+1) = input.performance.eliminate_DL.default;
        end
    end
    
    
    assembly_input.performance.DL_singularity_removal = logical([]);
    for field = string({Mesh.name})
        if isfield(input.performance.DL_singularity_removal,field)
            assembly_input.performance.DL_singularity_removal(end+1) = input.performance.DL_singularity_removal.(field);
        else
            assembly_input.performance.DL_singularity_removal(end+1) = input.performance.DL_singularity_removal.default;
        end
    end
    
    %     assembly_input.performance.DL_singularity_removal = input.performance.DL_singularity_removal;
    assembly_input.rotating_flagellum = input.rotating_flagellum;
    assembly_input.Tail.motorBC = input.Tail.motorBC;
    
    assembly_input.accuracy.ignore_interaction = input.accuracy.ignore_interaction;
    assembly_input.accuracy.integration_tol.stokeslet = input.accuracy.integration_tol.stokeslet;
    assembly_input.accuracy.integration_tol.stresslet = input.accuracy.integration_tol.stresslet;
    assembly_input.accuracy.integration_tol.force = input.accuracy.integration_tol.force;
    assembly_input.accuracy.integration_tol.torque = input.accuracy.integration_tol.torque;
    % currently the order of dino submeshes when loaded is not
    % carefully controlled so doing this here makes sure vector of eps2
    % matches order of Mesh.name
    assembly_input.accuracy.eps2 = [];
    for field = string({Mesh.name})
        if isfield(input.accuracy.epsilon,field)
            assembly_input.accuracy.eps2(end+1) = input.accuracy.epsilon.(field)^2;
        else
            assembly_input.accuracy.eps2(end+1) = input.accuracy.epsilon.default^2;
        end
    end
    
    if ~any([Mesh.is_mesh_rigid])
        assembly_input.performance.rigid_body_matrix_rotation = false; % regardless of what is requested in input, there's no where to use this shortcut
    else
        assembly_input.performance.rigid_body_matrix_rotation = input.performance.rigid_body_matrix_rotation;
    end
    
    %     assembly_input.constants.multfactor = input.constants.multfactor;
    assembly_input.constants.mu = input.constants.mu;  %not actually needed for matrix assembly, but needed for yprime() during timestepping interpolation
    assembly_input.performance.debug_mode = input.performance.debug_mode;
    assembly_input.performance.numels_max = input.performance.numels_max;
    assembly_input.performance.verbose = input.performance.verbose;
    
    
    
    switch input.problemtype
        case "resistance"
            assembly_input.Tail.motor_torque = NaN;
            
            flowcases = {'x','y','z','rx','ry','rz'}; %translations and rotations for each axis
            
            %             assembly_input.skip_rigid_integrals = false;
            
            %forced rotations will occur around refpoint
            % should update this to always use same refpoint, possibly the origin, to avoid confusion as to how D.center (center of diffusion) is defined!
            %         if strcmp(input.bugtype,'bacteria') && length(Mesh) == 2
            %             refpoint = Mesh(2).refpoints(:,1); %make sure refpoint is somewhere along the axis of motor rotation (actually this isn't necessary here....)
            %         else
            %             refpoint = Mesh(1).refpoints(:,1); % doesn't really matter where it is, so use body or whatever first submesh is
            %         end
            %
            
            for i = 1:length(flowcases)
                flowcase = flowcases{i};
                
                
                U0 = 50;  %microns / s    arbitrary, for calculation of friction coeffs
                Omega0 = 1;  %rad / s     arbitrary, for calculation of rotational friction coeffs
                
                BCs.resistance(i) = set_resistance_BCs(U0, Omega0, flowcase); %defines BCs struct with BCs(i).resistance.U and BCs(i).resistance.Omega vectors
                
                %                 [RHS(:,i)] = assemble_RHS(input,Mesh, matrix_props, BCs(i));
                
            end
            
            
            % need to add velocities to Mesh here in case we are including DL in matrix
            % assembly
            Mesh = compute_BCs_mex(BCs.resistance,Mesh,assembly_input); % adds field u = fixed frame velocities at each vert to each submesh
            % each Mesh(i).u is Nx3xM where M is length(flowcases), N is #verts
            [matrix_props] = gen_matrix_props(input,Mesh, node_parameters);
            
            
            if matrix_props.n_unknown_u > 0 && size(Mesh(1).u,3) > 1
                disp('multiple flowcases specified for free slip BC'); % not even sure why this is a problem but I seem to have thought so at one point
%                 pause
            end
            
            if recycle_A %don't need to run matrix_assembly_mex at all, we just ran mobility before and calculated everything we need here
                
                A = A_recycle(1:matrix_props.n_rows,1:matrix_props.n_cols); %copy traction integrals
                clear A_recycle
                A_force = A_force_recycle; %was output from prev mobility run
                A_torque = A_torque_recycle; %was output from prev mobility run
                
            else %need to assemble A from scratch, but can save results for reuse in mobility run
                %                 assembly_input.skip_traction_integrals = false;
                %                 assembly_input.bugtype = input.bugtype;
                %  [A,~,A_force,A_torque,~] = matrix_assembly_mexed(Mesh,matrix_props,assembly_input,skip_rigid_integrals); %A_motor_torque doesn't mean anything here, and don't need debug_info
                if input.performance.verbose
                    disp('Assembling matrix for resistance case');
                end
                resistance_assembly_tic = tic;
                
                % [A,~,A_force,A_torque,~] = matrix_assembly_mexed_orig(Mesh,matrix_props,assembly_input); %A_motor_torque doesn't mean anything here, and don't need debug_info
                %                 [ A, ~,A_force,A_torque,debug_info] = matrix_assembly_mex_wrapper( Mesh, matrix_props,assembly_input );
                
               
                [ A, A_force, A_torque, RHS, A_motor_torque] = matrix_assembly_mex_wrapper(Mesh,matrix_props,index_mapping,node_parameters,assembly_input);
                %  toc(resistance_assembly_tic);    % 1079 s for original matrix_assembly_mexed
                timings.resistance_matrix_assembly = toc(resistance_assembly_tic);
                if input.performance.verbose
                    disp(['matrix_assembly_mex_wrapper took ',num2str(timings.resistance_matrix_assembly)]);
                end
                
                A_recycle = A; %save for possible use in mobility case after this
                
                if input.performance.debug_mode
                    save([input.paths.dumpfolder,input.paths.namebase.full,'_debug_info','.mat'],'debug_info');
                end
                
            end
            
            
            %             clear fcoeffs solutions RHS
            %             resistance_solves_tic = tic;
            %             for i = 1:length(flowcases)
            %                 flowcase = flowcases{i};
            %
            %
            %                 U0 = 50;  %microns / s    arbitrary, for calculation of friction coeffs
            %                 Omega0 = 1;  %rad / s     arbitrary, for calculation of rotational friction coeffs
            %
            %                 BCs(i).resistance = set_resistance_BCs(U0, Omega0, flowcase); %defines BCs struct with U and Omega vectors
            %
            %                 [RHS(:,i)] = assemble_RHS(input,Mesh, matrix_props, BCs(i));
            %
            %             end
            
            
            [F] = matrix_solve(input, matrix_props, A, RHS); %solves 6 matrix equations and returns 6 solution vectors (in this case just traction) as columns of F
            
            for i = 1:length(flowcases)
                flowcase = flowcases{i};
                
                [forces] = integrate_traction(A_force,A_torque,F(:,i));  %integrates traction and returns forces and torques on entire object
                [fcoeffs_temp] = friction_coeffs(forces,BCs.resistance(i));  %computes translational and rotational friction coeffs
                
                solutions.(flowcase).f = F(:,i);
                solutions.(flowcase).forces = forces;
                switch flowcase
                    case {'x','y','z'}
                        solutions.(flowcase).U0 = U0;
                    case {'rx','ry','rz'}
                        solutions.(flowcase).Omega0 = Omega0;
                end
                
                
                
                switch flowcase
                    case 'x'
                        fcoeffs.translation(1) = fcoeffs_temp.translation(1);
                    case 'y'
                        fcoeffs.translation(2) = fcoeffs_temp.translation(2);
                    case 'z'
                        fcoeffs.translation(3) = fcoeffs_temp.translation(3);
                    case 'rx'
                        fcoeffs.rotation(1) = fcoeffs_temp.rotation(1);
                    case 'ry'
                        fcoeffs.rotation(2) = fcoeffs_temp.rotation(2);
                    case 'rz'
                        fcoeffs.rotation(3) = fcoeffs_temp.rotation(3);
                end
                
            end %flowcases
            solutions.fcoeffs = fcoeffs;
            
            timings.resistance_matrix_solves = toc(resistance_solves_tic);
            
            clear A F forces fcoeffs_temp debug_info %A is by far the largest variable and perhaps not very useful later on
            %   save([dumpfolder,namebase,'_',flowcase,'.mat']);
            % save([dumpfolder,namebase,'.mat']);
            %   save([dumpfolder,namebase,prefix,'_','eps','_',num2str(epsilon),'.mat']);
            
            [D] = diffusivities(solutions);
            
            solutions_temp = solutions;  %for possible use after timestepping to get effective diffusivity
            for i = 1:length(flowcases)
                flowcase = flowcases{i};
                solutions_temp.(flowcase) = rmfield(solutions_temp.(flowcase),'f');  %don't need this for effective diffusivity and it takes up space
            end
            
        case 'mobility'
            switch input.bugtype
                case 'bacteria'
                    
                    switch input.Tail.motorBC
                        case 'torque'
                            BCs.mobility.motor_torque = input.Tail.motor_torque;
                            %                                 y0 = [Mesh(1).refpoints(:,1); [0 0 0 0]' ]; %use body refpoint by convention
                            y0 = zeros(7,1); % with fixed timestepping algorithm, y(1:3) represents a translation, which starts at zero
                        case 'freq'
                            BCs.mobility.motor_freq = input.tail.motor_freq;
                            %                                 y0 = [Mesh(1).refpoints(:,1); [0 0 0 ]' ];
                            y0 = zeros(6,1); % with fixed timestepping algorithm, y(1:3) represents a translation, which starts at zero
                    end
                case 'dino'
                    
                    y0 = zeros(6,1); % with fixed timestepping algorithm, y(1:3) represents a translation, which starts at zero
            end
            
            
            
            
            if input.performance.kinematics_interpolation
                interp_y_tic = tic;
                %                     try
                interp_y_test2;
                %                 interp_y;
                %                     catch
                %                         disp('Error in interp_y.  Continuing to next sweep iter.');
                %                         close
                %                         continue
                %                     end
                timings.interpolation = toc(interp_y_tic);
                %% unroll all fields and vector components of last and most accurate interpolant
                cc = 0;  clear best_interpolant
                fields = fieldnames(interpolant(end));
                for f = 1:length(fields)  %each vector variable, e.g. U, Omega, omega
                    for i = 1:length( interpolant(end).(fields{f}) )  %components of each variable, e.g U(1:3), Omega(1:3)
                        cc = cc+1;
                        best_interpolant(cc) = interpolant(end).(fields{f})(i);
                    end
                end
                
            else %can't do kinematics interpolation.  instead, must redo matrix assembly and matrix solve for every timestep.
                % currently not updated / broken!!!
                
                %                 global n_out solf
                %                 n_out = 1; %index for solf output
                %
                %                 solf = NaN(50,n_col*3+1+6+1+1);  %also time rotationrate tailtorque
                %
                %
                %                 skip_rigid_integrals = false;  %do full solve just once
                %                 [A0,Integral,Abserror,Numevals,Exitflag,coeffs0] = matrix_assembly_adaptive_vectorized_free_mexed(Mesh,n_rows,n_cols,n_col,Col_inds,eps2,progmon,nthreads,integration_tol,debug_output,'mobility',parameters(c).tail.motorBC,parameters(c).constants.multfactor,skip_rigid_integrals,ignore_interaction);
                %
                %
                %                 fun = @(t,y) yprime(t,y,Mesh,n_rows,n_cols,n_col,Col_inds,eps2,progmon,nthreads,integration_tol,parameters(c),A0,coeffs0,timerval,[0 T],true,ignore_interaction);
                
            end
            
            clear A0 A
            
            
            
    end %problemtype
    
    timings.total = toc(total_time_tic);
    if input.performance.verbose
        disp(['Sweep iter ',num2str(sweep_i),' of ',num2str(length(Inputs)),' took ',num2str(timings.total)]);
    end
    clear temp* x y
    temp = who;
    temp = setdiff(temp,{'A_recycle','A_force_recycle','A_torque_recycle', 'Inputs', 'A0', 'debug_info','Results','temp'}); %don't want to save A_recycle in the .mat file
    allvars = whos;
    tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
    tosave = {allvars(tosave).name};
    temp = intersect(temp,tosave);
    
    return
    save([input.paths.dumpfolder,input.paths.namebase.full,'_dump','.mat'],temp{:});
    
    
    clear solutions Solutions
    
    
    
    if strcmp(input.problemtype,'mobility')  && input.do_timestepping
        
        
        %  y0 = y(1:6);  input.phase_speed = Metadata_temp.geom.phase_speed;
        clear dump
        dump.y0 = y0;  dump.best_interpolant = best_interpolant;
        if input.performance.verbose
            disp('Beginning timestepping.');
        end
        timestepping_solution0 = [];  refpoint0 = Mesh(1).refpoints(:,1);
        [fits, timestepping_solution, timings] = timestepping(input,dump,timestepping_solution0,refpoint0,Mesh);  %timestepping dump file is saved internally, contains timestepping_solution and fits
        %             if input.performance.verbose
        %                 disp(['Timestepping took ',num2str(timings.timestepping / 60), ' min']);
        %             end
        
        
        %%
        
        if strcmp(input.bugtype,'bacteria')
            interpolant = interpolant(end).omega;  %last and most accurate interpolant for tail rotation rate omega
            
            [avg_omega] = compute_avg_omega(interpolant);
            
            %in addition to storing in memory for later inclusion into aggregate results file, immediately save avg_omega and motor_torque inside timestepping dump file
            m = matfile([input.paths.dumpfolder,input.paths.namebase.full,'_timestepping','.mat'],'Writable',true);
            m.avg_omega = avg_omega;
            m.motor_torque = input.tail.motor_torque;
            
            fits.avg_swimming_axis = calc_avg_swimming_axis(fits, timestepping_solution, Mesh); %outputs direction of avg swimming direction in body frame
            
            m.fits = fits;
            
            
            
        end %if bacteria
        
    end  %if mobility, do timestepping right now
    
    if strcmp(input.bugtype,'bacteria')
        
        save_Results;
    end
    
    
end  %parameter space sweep

% toc(inputs_loop_tic)

clear A_recycle A A0 debug_info