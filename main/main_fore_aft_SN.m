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
    
    clear Mesh Metadata timings BCs
    
    %% all this is too ensure we can reuse part of A matrix between forced and freeswim, in conjunction with randomize_verts
    if sweep_i > 1  %can't reuse anything if this is the first sweep iteration, duh
        [diff_fields] = struct_diff(Inputs(sweep_i-1),Inputs(sweep_i));  %if we can reuse, prev inputs should be identical except for a few things:
        input_diffs = setdiff(diff_fields,{'full','problemtype'});
    else
        input_diffs = NaN;  %anything that's not empty works
    end
    
    if isempty(input_diffs)  %we just did forced and are now doing freeswim or vice versa
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
            if strcmp(input.problemtype,'freeswim') || ( strcmp(input.problemtype,'forced') && input.include_tail )
                submeshes = {'body','tail'};  %which submeshes to load  (always list body first by convention)
            else
                submeshes = {'body'};
            end
        case 'dino'
            % temp = {'Body','Transverse','Tail','Wingtip'};
            %             submeshes = {'body','transverse','tail'};
            %  submeshes = temp(logical(input.potatohead));  %only load submeshes that will be simulated - not sure if this will break shat for old cases.....
            
            % changed to automatically figure out which submeshes were used
            % during compute BCs internally, and load all of them
            [Mesh_files] = get_dino_mesh_files(input.paths.datfolder);  %Mesh_files is a list of the mesh and Metadata filenames for all the phase angles
        case 'sheet'
            submeshes = {'sheet'};
            [Mesh_files] = get_sheet_mesh_files(input.paths.datfolder);
    end
    

    %for monoflagellated bacteria model, Mesh(1) is the body, Mesh(2) is the tail
    if strcmp(input.bugtype,'bacteria')   || (  strcmp(input.bugtype,'dino') && strcmp(input.problemtype,'forced')  )  %dino freeswim case has meshes loaded again later anyway in interp_y
        bad_sweep_i(sweep_i) = false;
        for si = 1:length(submeshes)
            meshname = [input.paths.datfolder,input.paths.namebase.(submeshes{si}),'.dat'];
            if input.performance.verbose
                disp(['Loading ',submeshes{si},' mesh     ',meshname]);
            end
            
            
            try
                [temp_mesh, Metadata(si), persistant_data(si).rand_inds] = load_mesh_wrapper(meshname, input, persistant_data(si).rand_inds, input_diffs);  %old rand_inds goes in, new rand_inds comes out (which might be the same as it was)
            catch
                [mesh_succeed] = gen_mesh_wrapper(input);
                if ~mesh_succeed
                    stopafra
                end
                [temp_mesh, Metadata(si), persistant_data(si).rand_inds] = load_mesh_wrapper(meshname, input, persistant_data(si).rand_inds, input_diffs);  %old rand_inds goes in, new rand_inds comes out (which might be the same as it was)
            end
            
            if ( isfield(Metadata(si).mesh,'meshing_succeeded') && ~Metadata(si).mesh.meshing_succeeded ) || isempty(temp_mesh)
                bad_sweep_i(sweep_i) = true;
                break %don't bother trying to load other submeshes
            end
            
            temp_mesh.name = submeshes{si};
            temp_mesh.orientation = [1 0 0; 0 1 0; 0 0 1;]';  %initial orientation vectors are along x, y, z axes
            
            Mesh(si) = temp_mesh;  clear temp_mesh
            
        end
        
        if bad_sweep_i(sweep_i)
            continue
        end
        
    end
    
    
    if strcmp(input.bugtype,'bacteria') && length(Mesh) >= 2  %make sure no elems or verts are mistaken as shared with body, since meshes were created separately unlike dino case
        Mesh(2).indices.orig.elem = Mesh(2).indices.orig.elem + max(Mesh(1).indices.orig.elem);
        Mesh(2).indices.orig.vert = Mesh(2).indices.orig.vert + max(Mesh(1).indices.orig.vert);
        Mesh(2).elems = Mesh(2).elems + max(Mesh(1).indices.orig.vert);  %still using orig vert labels here, until renumber_Mesh below
    elseif (  strcmp(input.bugtype,'dino') && strcmp(input.problemtype,'forced')  )  % tail isn't connected to body/transverse so could have overlapping indices
        stoperoo
        % this prolly doesn't work anymore and needs to be updated....
        
        Mesh(3).indices.orig.elem = Mesh(3).indices.orig.elem + max([Mesh(1).indices.orig.elem; Mesh(2).indices.orig.elem]);
        Mesh(3).indices.orig.vert = Mesh(3).indices.orig.vert + max([Mesh(1).indices.orig.vert; Mesh(2).indices.orig.vert]);
        Mesh(3).elems = Mesh(3).elems + max([Mesh(1).indices.orig.vert; Mesh(2).indices.orig.vert]);  %still using orig vert labels here, until renumber_Mesh below
    end
    
    if strcmp(input.bugtype,'bacteria')  || (  strcmp(input.bugtype,'dino') && strcmp(input.problemtype,'forced')  )
        [Mesh] = global_inds(Mesh);
        [Mesh] = renumber_Mesh(Mesh);
        % create smaller input struct with just necessary parameters to reduce mex
        % pain
        temp_input.performance.nthreads = input.performance.nthreads;
        temp_input.accuracy.integration_tol.area = input.accuracy.integration_tol.area;
        temp_input.accuracy.integration_tol.centroid = input.accuracy.integration_tol.centroid;
        temp_input.accuracy.integration_tol.volume = input.accuracy.integration_tol.volume;
        [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
    end
    
    if strcmp(input.bugtype,'bacteria')
        
        if strcmp(input.body.shape,'curved_rod')
            
            if input.body.AR(1) == 1 %sphere degenerate case
                Mesh(1).refpoints = [0 0 0]';   %center
            elseif input.body.AR(2) == 0  %not a sphere, but a straight rod degenerate case
                Mesh(1).refpoints = [0 0 0]';   %center
            else  %actually a general curved rod
                Mesh(1).refpoints =  [2*Metadata(1).geom.radius_curv*sin(Metadata(1).geom.nturns / 2 *pi) * cos(Metadata(1).geom.nturns / 2 *pi); ...
                    2*Metadata(1).geom.radius_curv*(sin(Metadata(1).geom.nturns / 2 *pi))^2;
                    0 ];   %center of centerline
            end
            
        else
            Mesh(1).refpoints = [0 0 0]';   %probably center
        end
        
        if length(Mesh) == 2 && input.include_tail
            Mesh(2).refpoints = [0 0 0]';  % will be used for torque calculations - OK as long as somewhere along motor axis
        end
    end
    
    
    if  strcmp(input.bugtype,'bacteria') && input.include_tail && input.tail.rotate_tail %to check for possible effect of tail phase angle on calculated diffusivities   %strcmp(input.problemtype,'forced') &&
        % do this even for freeswim so that if the tail was rotated for
        % forced, it will get rotated for freeswim and thus part of A
        % matrix can be reused.  Not that we're likely to ever do
        % freeswimming runs in conjunction with this option anyway....
        [Mesh(2)] = rotateMesh(Mesh(2), [ input.tail.rotate_tail_angle  0  0]' );  %rotate tail around x axis
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
            switch input.body.shape
                case 'ellipsoid' %origin for body located at body center
                    Mesh(2) = shiftMesh(Mesh(2),[-2*input.tail.radius - Metadata(1).geom.a 0 0]);
                case 'capsule' %origin for body located at body center
                    Mesh(2) = shiftMesh(Mesh(2),[-2*input.tail.radius - (Metadata(1).geom.caps_c + Metadata(1).geom.height/2) 0 0]);
                case 'curved_rod'  %origin for body located at center of spherical cap
                    if input.body.AR(1) == 1 %sphere degenerate case
                        Mesh(2) = shiftMesh(Mesh(2),[-2*input.tail.radius - Metadata(1).geom.radius 0 0]);
                    elseif input.body.AR(2) == 0  %not a sphere, but a straight rod degenerate case
                        Mesh(2) = shiftMesh(Mesh(2),[-2*input.tail.radius - (Metadata(1).geom.radius + Metadata(1).geom.height/2) 0 0]);
                    else  %actually a general curved rod
%                         if ~exist('tail_shift','var')
%                             tail_shift = 0;
%                         end
                        
                        Mesh(2) = shiftMesh(tail_mesh,[-2*input.tail.radius - (Metadata(1).geom.radius) - tail_shift 0 0]);
                    end
                    
            end
        else
            break_switch = true;
            
        end
        
        if ( strcmp(input.bugtype,'bacteria') && length(Mesh) > 1 )  || (  strcmp(input.bugtype,'dino') && strcmp(input.problemtype,'forced')  )
            intersection_check_tic = tic;
            if input.body.AR(2) <= 0.85
                disp('skipping intersections check');
                is_intersected = false;
            else
                nrefines = [2 1];
                [is_intersected] = submesh_intersections(Mesh,input.accuracy.check_intersections_tolerance, input.accuracy.check_intersections_n_angles,false,input.performance.nthreads,nrefines);  % checks for self-intersection at each angle in parallel
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
        if (strcmp(input.problemtype,'forced') && input.include_tail)
            Centroid_overall = sum(repmat([Mesh.Volume],3,1) .* [Mesh.Centroid],2) ./ sum([Mesh.Volume]); %Centroid of entire mesh
            Mesh = shiftMesh(Mesh, - Centroid_overall); %after this point, motor axis no longer coincides with x-axes, so don't try to rotateMesh the tail!
        end
        
        
        
        if strcmp(input.bugtype,'bacteria')  || (  strcmp(input.bugtype,'dino') && strcmp(input.problemtype,'forced')  )
            clear matrix_props
            % matrix_props.n_col = sum([Mesh.n_vert]);  %number collocation points in A matrix, could be greater than number of points at which to determine traction, in theory
            %matrix_props.n_col = Mesh(end).indices.glob.vert( find( ~isnan(Mesh(end).indices.glob.vert), 1, 'last' ) );  %highest global index must be in last submesh, last nonNaN entry since it must monotonically increase
            matrix_props.n_col = Mesh(end).indices.glob.unq_bounds.vert(2);
            
            
            if strcmp(input.problemtype,'freeswim')
                switch input.bugtype
                    case 'bacteria'
                        switch input.tail.motorBC
                            case 'freq'
                                matrix_props.n_rows = matrix_props.n_col * 3 + 6;  %usual unknowns plus 3 translation components and 3 rotation components of body
                            case 'torque'
                                matrix_props.n_rows = matrix_props.n_col * 3 + 7;  %usual unknowns plus 3 translation components and 3 rotation components of body + rotation rate of tail
                        end
                    case 'dino'
                        matrix_props.n_rows = matrix_props.n_col * 3 + 6;  %usual unknowns plus 3 translation components and 3 rotation components of body
                end
            else
                matrix_props.n_rows = matrix_props.n_col * 3;
            end
            matrix_props.n_cols = matrix_props.n_rows; %number of columns in A matrix
            
            matrix_props.Col_inds = save_Col_inds_mexed(Mesh,input.performance.nthreads);
            
        end  %if bacteria or forced dino
        
        
        % make smaller version of input struct with just variables needed to assemble A
        % matrix, to reduce pain of re-mexing if input struct is modified
        clear assembly_input
        assembly_input.performance.nthreads = input.performance.nthreads;
        assembly_input.problemtype = input.problemtype;
        assembly_input.tail.motorBC = input.tail.motorBC;
        assembly_input.ignore_interaction = input.ignore_interaction;
        assembly_input.accuracy.integration_tol.traction = input.accuracy.integration_tol.traction;
        assembly_input.accuracy.integration_tol.force = input.accuracy.integration_tol.force;
        assembly_input.accuracy.integration_tol.torque = input.accuracy.integration_tol.torque;
        assembly_input.accuracy.eps2 = input.accuracy.eps2;
        assembly_input.constants.multfactor = input.constants.multfactor;
        assembly_input.constants.mu = input.constants.mu;  %not actually needed for matrix assembly, but needed for yprime() during timestepping interpolation
        assembly_input.performance.debug_mode = input.performance.debug_mode;
        assembly_input.performance.numels_max = input.performance.numels_max;
        assembly_input.performance.verbose = input.performance.verbose;
        
     
        
        switch input.problemtype
            case 'forced'
                flowcases = {'x','y','z','rx','ry','rz'}; %translations and rotations for each axis
                
                assembly_input.skip_rigid_integrals = false;
                
                %             assembly_input.skip_rigid_integrals = true;
                %             assembly_input.problemtype = 'freeswim';
                
                if recycle_A %don't need to run matrix_assembly_mex at all, we just ran freeswim before and calculated everything we need here
                    
                    A = A_recycle(1:matrix_props.n_rows,1:matrix_props.n_cols); %copy traction integrals
                    clear A_recycle
                    A_force = A_force_recycle; %was output from prev freeswim run
                    A_torque = A_torque_recycle; %was output from prev freeswim run
                    
                else %need to assemble A from scratch, but can save results for reuse in freeswim run
                    assembly_input.skip_traction_integrals = false;
                    assembly_input.bugtype = input.bugtype;
                    %  [A,~,A_force,A_torque,~] = matrix_assembly_mexed(Mesh,matrix_props,assembly_input,skip_rigid_integrals); %A_motor_torque doesn't mean anything here, and don't need debug_info
                    if input.performance.verbose
                        disp('Assembling matrix for forced case');
                    end
                    forced_assembly_tic = tic;
                    
                    % [A,~,A_force,A_torque,~] = matrix_assembly_mexed_orig(Mesh,matrix_props,assembly_input); %A_motor_torque doesn't mean anything here, and don't need debug_info
                    [ A, ~,A_force,A_torque,debug_info] = matrix_assembly_mex_wrapper( Mesh,matrix_props,assembly_input );
                    
                  %  toc(forced_assembly_tic);    % 1079 s for original matrix_assembly_mexed
                    timings.forced_matrix_assembly = toc(forced_assembly_tic);
                    if input.performance.verbose
                        disp(['matrix_assembly_mex_wrapper took ',num2str(timings.forced_matrix_assembly)]);
                    end
                    
                    A_recycle = A; %save for possible use in freeswim case after this
                    
                    if input.performance.debug_mode
                        save([input.paths.dumpfolder,input.paths.namebase.full,'_debug_info','.mat'],'debug_info');
                    end
                    
                end
                
                
                clear fcoeffs solutions
                forced_solves_tic = tic;
                for i = 1:length(flowcases)
                    flowcase = flowcases{i};
                    
                    
                    U0 = 50;  %microns / s    arbitrary, for calculation of friction coeffs
                    Omega0 = 1;  %rad / s     arbitrary, for calculation of rotational friction coeffs
                    
                    BCs.forced = set_forced_BCs(U0, Omega0,flowcase); %defines BCs struct with U and Omega vectors
                    
                    [RHS] = assemble_RHS(input,Mesh, matrix_props, BCs);
                    [f] = matrix_solve(input, matrix_props, A, RHS); %solves matrix equation and returns solution vector, in this case just traction
                    [forces] = integrate_traction(A_force,A_torque,f);  %integrates traction and returns forces and torques on entire object
                    [fcoeffs_temp] = friction_coeffs(forces,BCs.forced);  %computes translational and rotational friction coeffs
                    
                    solutions.(flowcase).f = f;
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
                
                timings.forced_matrix_solves = toc(forced_solves_tic);
                
                clear A f forces fcoeffs_temp debug_info %A is by far the largest variable and perhaps not very useful later on
                %   save([dumpfolder,namebase,'_',flowcase,'.mat']);
                % save([dumpfolder,namebase,'.mat']);
                %   save([dumpfolder,namebase,prefix,'_','eps','_',num2str(epsilon),'.mat']);
                
                [D] = diffusivities(solutions);
                
                solutions_temp = solutions;  %for possible use after timestepping to get effective diffusivity
                for i = 1:length(flowcases)
                    flowcase = flowcases{i};
                    solutions_temp.(flowcase) = rmfield(solutions_temp.(flowcase),'f');  %don't need this for effective diffusivity and it takes up space
                end
                
            case 'freeswim'
                switch input.bugtype
                    case 'bacteria'
                        switch input.tail.motorBC
                            case 'torque'
                                BCs.freeswim.motor_torque = input.tail.motor_torque;
                                y0 = [Mesh(1).refpoints(:,1); [0 0 0 0]' ]; %use body refpoint by convention
                            case 'freq'
                                BCs.freeswim.motor_freq = input.tail.motor_freq;
                                y0 = [Mesh(1).refpoints(:,1); [0 0 0 ]' ];
                        end
                    case 'dino'
                        % everything is done in interp_y, theoretically
                end
                
                
                
                
                if input.performance.kinematics_interpolation
                    interp_y_tic = tic;
                    try
                        interp_y;
                    catch
                        disp('Error in interp_y.  Continuing to next sweep iter.');
                        close
                        continue
                    end
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
                    %                 [A0,Integral,Abserror,Numevals,Exitflag,coeffs0] = matrix_assembly_adaptive_vectorized_free_mexed(Mesh,n_rows,n_cols,n_col,Col_inds,eps2,progmon,nthreads,integration_tol,debug_output,'freeswim',parameters(c).tail.motorBC,parameters(c).constants.multfactor,skip_rigid_integrals,ignore_interaction);
                    %
                    %
                    %                 fun = @(t,y) yprime(t,y,Mesh,n_rows,n_cols,n_col,Col_inds,eps2,progmon,nthreads,integration_tol,parameters(c),A0,coeffs0,timerval,[0 T],true,ignore_interaction);
                    
                end
                
                clear A0 A
                
                try
                    saveas(input.output.interpolation.fignum,[input.paths.dumpfolder,input.paths.namebase.full,'_interpolation','.fig']);
                catch
                end
                
        end %problemtype
        
        timings.total = toc(total_time_tic);
        if input.performance.verbose
            disp(['Sweep iter ',num2str(sweep_i),' of ',num2str(length(Inputs)),' took ',num2str(timings.total)]);
        end
        clear temp* x y
        temp = who;
        temp = setdiff(temp,{'A_recycle','A_force_recycle','A_torque_recycle', 'Inputs', 'A0', 'debug_info','Results','temp'}); %don't want to save A_recycle in the .mat file
        
        save([input.paths.dumpfolder,input.paths.namebase.full,'_dump','.mat'],temp{:});
        
        
        clear solutions Solutions
        
        
        
        if strcmp(input.problemtype,'freeswim')  && input.do_timestepping
            
            
            %  y0 = y(1:6);  input.phase_speed = Metadata_temp.geom.phase_speed;
            
            dump.y0 = y0;  dump.best_interpolant = best_interpolant;
            if input.performance.verbose
                disp('Beginning timestepping.');
            end
            timestepping_solution0 = [];
            [fits, timestepping_solution, timings] = timestepping(input,dump,timestepping_solution0);  %timestepping dump file is saved internally, contains timestepping_solution and fits
            if input.performance.verbose
                disp(['Timestepping took ',num2str(timings.timestepping / 60), ' min']);
            end
            
            if strcmp(input.bugtype,'bacteria')
                interpolant = interpolant(end).omega;  %last and most accurate interpolant for tail rotation rate omega
                
                [avg_omega] = compute_avg_omega(interpolant);
                
                %in addition to storing in memory for later inclusion into aggregate results file, immediately save avg_omega and motor_torque inside timestepping dump file
                m = matfile([input.paths.dumpfolder,input.paths.namebase.full,'_timestepping','.mat'],'Writable',true);
                m.avg_omega = avg_omega;
                m.motor_torque = input.tail.motor_torque;
                
                
                
                
            end %if bacteria
            
        end  %if freeswim, do timestepping right now
        
        if strcmp(input.bugtype,'bacteria')
            
            save_Results_fore_aft_SN;
        end
        
        
    end  %parameter space sweep
    
    % toc(inputs_loop_tic)
    
    clear A_recycle A A0 debug_info