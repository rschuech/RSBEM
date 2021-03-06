

%clearvars -except Inputs
inputs_loop_tic = tic;
%if run from Condor, current directory will be temporary folder where
%Condor has placed wrapper.m
%cd E:/Hull/code/;
%cd C:\RD\bem\code\;
addpath(genpath('./')); %add all subfolders to path so that subfunctions will be found


%bad_sweep_i = false(1,length(Inputs));  %keep track of which runs were bad due to not having mesh files
for sweep_i = 1:length(Inputs)
    disp(['On sweep iter ',num2str(sweep_i),' of ',num2str(length(Inputs)),' (not including parallel timestepping)']);
    total_time_tic = tic;
    input = Inputs(sweep_i);
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
    
    if ~exist('persistant_data','var')
        persistant_data(1).rand_inds = [];  %use placeholder to start; this will hold rand_inds for each submesh and maybe other things to save between sweep iters later
        persistant_data(2).rand_inds = [];
    end
    %%
    
    %% load body mesh
    
    
    switch input.bugtype
        case 'bacteria'
            if strcmp(input.problemtype,'freeswim') || ( strcmp(input.problemtype,'forced') && input.include_tail )
                submeshes = {'body','tail'};  %which submeshes to load  (always list body first by convention)
            else
                submeshes = {'body'};
            end
        case 'dino'
            submeshes = {'body','transverse','tail'};
    end
    %for monoflagellated bacteria model, Mesh(1) is the body, Mesh(2) is the tail
      bad_sweep_i(sweep_i) = false;
    for si = 1:length(submeshes)
        meshname = [input.paths.datfolder,input.paths.namebase.(submeshes{si}),'.dat'];
        
        disp(['Loading ',submeshes{si},' mesh     ',meshname]);
        [temp_mesh, Metadata(si), persistant_data(si).rand_inds] = load_mesh_wrapper(meshname,input,persistant_data(si).rand_inds,input_diffs);  %old rand_inds goes in, new rand_inds comes out (which might be the same as it was)
        
        if ( isfield(Metadata(si).mesh,'meshing_succeeded') && ~Metadata(si).mesh.meshing_succeeded ) || isempty(temp_mesh)
            bad_sweep_i(sweep_i) = true;
            break
        end
        
        temp_mesh.name = submeshes{si};
        temp_mesh.orientation = [1 0 0; 0 1 0; 0 0 1;]';  %initial orientation vectors are along x, y, z axes
        
        
        
        Mesh(si) = temp_mesh;  clear temp_mesh
        
    end
    
     if bad_sweep_i(sweep_i)
         continue
     end
    
    if strcmp(input.bugtype,'bacteria')  %make sure no elems or verts are mistaken as shared with body, since meshes were created separately unlike dino case
        Mesh(2).indices.orig.elem = Mesh(2).indices.orig.elem + max(Mesh(1).indices.orig.elem);
         Mesh(2).indices.orig.vert = Mesh(2).indices.orig.vert + max(Mesh(1).indices.orig.vert);
         Mesh(2).elems = Mesh(2).elems + max(Mesh(1).indices.orig.vert);  %still using orig vert labels here, until renumber_Mesh below
    end
    
    
    [Mesh] = global_inds(Mesh);
    [Mesh] = renumber_Mesh(Mesh);
    % create smaller input struct with just necessary parameters to reduce mex
    % pain
    temp_input.performance.nthreads = input.performance.nthreads;
    temp_input.accuracy.integration_tol.area = input.accuracy.integration_tol.area;
    temp_input.accuracy.integration_tol.centroid = input.accuracy.integration_tol.centroid;
    temp_input.accuracy.integration_tol.volume = input.accuracy.integration_tol.volume;
    [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
    
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
        
         Mesh(2).refpoints = [0 0 0]';  % will be used for torque calculations - OK as long as somewhere along motor axis

    end
    
    
    if  input.tail.rotate_tail %to check for possible effect of tail phase angle on calculated diffusivities   %strcmp(input.problemtype,'forced') &&
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
    
  
    
    %%
    if strcmp(input.bugtype,'bacteria')
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
                    Mesh(2) = shiftMesh(Mesh(2),[-2*input.tail.radius - (Metadata(1).geom.radius) 0 0]);
                end
        end
        
        
        intersection_check_tic = tic;
        [is_intersected] = submesh_intersections(Mesh,input.accuracy.check_intersections_tolerance, input.accuracy.check_intersections_n_angles,false,input.performance.nthreads);  % checks for self-intersection at each angle in parallel
        timings.intersection_check = toc(intersection_check_tic);
        
        if is_intersected
             bad_sweep_i(sweep_i) = true;
            load(input.paths.intersections_file,'intersections'); %contains one variable called intersections, a list of non-working input structs
            if isempty(intersections)
                intersections = input;
            else
                intersections(end+1) = input;  %should contain all the info we'll ever need on the current geometry combination that doesn't work
            end
            save(input.paths.intersections_file,'intersections');
            disp('Tail intersects Body; aborting.');
            continue %with Inputs parameter sweep loop
        else
            disp('Intersections test passed.');
            continue
        end
        
        %if computing friction coeffs, shift whole thing so that fixed ref
        %frame origin is at centroid (though this isn't really necessary)
        if (strcmp(input.problemtype,'forced') && input.include_tail)
            Centroid_overall = sum(repmat([Mesh.Volume],3,1) .* [Mesh.Centroid],2) ./ sum([Mesh.Volume]); %Centroid of entire mesh
            Mesh = shiftMesh(Mesh, - Centroid_overall); %after this point, motor axis no longer coincides with x-axes, so don't try to rotateMesh the tail!
        end
        
    end %if bugtype == bacteria
    
   
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
  
    
    
    
    
    switch input.problemtype
        case 'forced'
            flowcases = {'x','y','z','rx','ry','rz'}; %translations and rotations for each axis
            
            assembly_input.skip_rigid_integrals = false;
            
            if recycle_A %don't need to run matrix_assembly_mex at all, we just ran freeswim before and calculated everything we need here
                
                A = A_recycle(1:matrix_props.n_rows,1:matrix_props.n_cols); %copy traction integrals
                A_force = A_force_recycle; %was output from prev freeswim run
                A_torque = A_torque_recycle; %was output from prev freeswim run
                
            else %need to assemble A from scratch, but can save results for reuse in freeswim run
                assembly_input.skip_traction_integrals = false;
                  assembly_input.bugtype = input.bugtype;
                %  [A,~,A_force,A_torque,~] = matrix_assembly_mexed(Mesh,matrix_props,assembly_input,skip_rigid_integrals); %A_motor_torque doesn't mean anything here, and don't need debug_info
                disp('Assembling matrix for forced case');
                forced_assembly_tic = tic;
                
                % [A,~,A_force,A_torque,~] = matrix_assembly_mexed_orig(Mesh,matrix_props,assembly_input); %A_motor_torque doesn't mean anything here, and don't need debug_info
                [ A, ~,A_force,A_torque,~] = matrix_assembly_mex_wrapper( Mesh,matrix_props,assembly_input );
                
                toc(forced_assembly_tic);    % 1079 s for original matrix_assembly_mexed
                timings.forced_matrix_assembly = toc(forced_assembly_tic);
                A_recycle = A; %save for possible use in freeswim case after this
                
            end
            
            
            clear fcoeffs
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
            timings.forced_matrix_solves = toc(forced_solves_tic);
            
            clear A f forces fcoeffs_temp %A is by far the largest variable and perhaps not very useful later on
            %   save([dumpfolder,namebase,'_',flowcase,'.mat']);
            % save([dumpfolder,namebase,'.mat']);
            %   save([dumpfolder,namebase,prefix,'_','eps','_',num2str(epsilon),'.mat']);
            
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
                    y0 = [Mesh(1).refpoints(:,1); [0 0 0 ]' ];
            end
            
            
            if input.performance.kinematics_interpolation
                interp_y_tic = tic;
                interp_y;
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
            
            saveas(input.output.interpolation.fignum,[input.paths.dumpfolder,input.paths.namebase.full,'_interpolation','.fig']);
            
    end %problemtype
    
    timings.total = toc(total_time_tic);
    
    temp = who;
    temp = setdiff(temp,{'A_recycle','A_force_recycle','A_torque_recycle', 'Inputs', 'A0'}); %don't want to save A_recycle in the .mat file
    
    save([input.paths.dumpfolder,input.paths.namebase.full,'_dump','.mat'],temp{:});
    clear solutions Solutions
    
end  %parameter space sweep

toc(inputs_loop_tic)
clear A_recycle
