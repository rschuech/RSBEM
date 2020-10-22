function [Ax_1,Ay_1,Az_1,  Ax_2,Ay_2,Az_2,  Ax_3,Ay_3,Az_3,  Ax_4,Ay_4,Az_4,    A_freeswim_rows,A_freeswim_cols,A_motor_torque,A_force,A_torque,  debug_info] = matrix_assembly_mex(Mesh,matrix_props,assembly_input)

%assembles all entries of A matrix for both forced and force-free cases
%for free-swim and a freq motorBC , A_motor_torque can be multiplied by f to give the motor
%torque
%for all cases but particularly the forced case, A_force and A_torque can
%be multiplied by f to give the total force and torque on the object
%if debug_mode is true, debug_info can be output

%only produces Ax, Ay, Az horizontal sections of full A to avoid mex error
%due to A being too big.  As a workaround, we insert Ax, Ay, Az into full A in
%matrix_assembly_mex_wrapper, in native Matlab, which has no dumb
%limitation on matrix size.  Both versions appear to take exactly the same
%amount of time.

if strcmp(assembly_input.bugtype,'bacteria') && length(Mesh) == 2
    motor_orientation = Mesh(2).orientation(:,1);  % updated to allow for tails attached at an angle - motor rotation axis should always be the centerline axis of the tail
    refpoint = Mesh(2).refpoints(:,1); %make sure refpoint is somewhere along the axis of motor rotation for free swimming problems involving a motor torque condition
else
    motor_orientation = [];  %appease Coder
    refpoint = Mesh(1).refpoints(:,1); % doesn't really matter where it is, so use body or whatever first submesh is
end


coder.extrinsic('tic');
coder.extrinsic('toc');
coder.extrinsic('num2str');

debug_info.abs_error = [];
debug_info.rel_error = [];
debug_info.fun_evals = [];
debug_info.flag = [];
debug_info.min_value = [];
debug_info.max_value = [];
% debug_info.time = [];


if length(Mesh) > 1 %multiple submeshes input, use global indices
    tot_elems = Mesh(end).indices.glob.unq_bounds.elem(2);
else %only one submesh input, use local indices
    tot_elems = Mesh(1).n_elem;
end

tic
%#codegen
%%
%preallocate A matrix and temporary submatrices needed for parallel assembly

%test anticipated size of Ax, Ay, Az against dumb coder max array size
%limit and if problemo, split Ax, Ay, Az further into sub-arrays
% splitting is based on columns to facilitate less parfor pain

%right now limited to a further division by 4 but this is easily expanded to
%whatever is needed

%numels_max = double(intmax('int32'));  %dumb Coder limit on numel of any array


if matrix_props.n_col * matrix_props.n_col*3 < assembly_input.performance.numels_max  %good enough to just use Ax, Ay, Az
    Ax_1 = zeros(matrix_props.n_col,matrix_props.n_col*3);  Ay_1 = Ax_1;  Az_1 = Ax_1;  %only includes tractions for all cases
    Ax_2 = [];  Ay_2 = [];  Az_2 = [];   Ax_3 = [];  Ay_3 = [];  Az_3 = [];   Ax_4 = [];  Ay_4 = [];  Az_4 = [];
    
elseif ceil(matrix_props.n_col * matrix_props.n_col*3 / 2)  < assembly_input.performance.numels_max  %what if we split Ax, Ay, Az into two, are they small enough then?
    Ax_1 = zeros(matrix_props.n_col,ceil(matrix_props.n_col*3 / 2));  Ay_1 = Ax_1;  Az_1 = Ax_1;  %only includes tractions for all cases
    Ax_2 = zeros(matrix_props.n_col,matrix_props.n_col*3 - ceil(matrix_props.n_col*3 / 2));  Ay_2 = Ax_2;  Az_2 = Ax_2;  %only includes tractions for all cases
    Ax_3 = [];  Ay_3 = [];  Az_3 = [];   Ax_4 = [];  Ay_4 = [];  Az_4 = [];
    if assembly_input.performance.verbose
        disp('Dividing A into 2 parts as Coder limitation workaround.')
    end
else %ceil(matrix_props.n_col * matrix_props.n_col*3 / 4)  < assembly_input.performance.numels_max  %what if we split Ax, Ay, Az into two, are they small enough then?
    Ax_1 = zeros(matrix_props.n_col,ceil(matrix_props.n_col*3 / 4));  Ay_1 = Ax_1;  Az_1 = Ax_1;  %only includes tractions for all cases
    Ax_2 = zeros(matrix_props.n_col,ceil(matrix_props.n_col*3 / 4));  Ay_2 = Ax_2;  Az_2 = Ax_2;  %only includes tractions for all cases
    Ax_3 = zeros(matrix_props.n_col,ceil(matrix_props.n_col*3 / 4));  Ay_3 = Ax_3;  Az_3 = Ax_3;  %only includes tractions for all cases
    Ax_4 = zeros( matrix_props.n_col ,matrix_props.n_col*3 - ceil(matrix_props.n_col*3 / 4)*3 );  Ay_4 = Ax_4;  Az_4 = Ax_4;  %only includes tractions for all cases
    if assembly_input.performance.verbose
        disp('Dividing A into 4 parts as Coder limitation workaround.')
    end
    %code will fail around here if we end up needing further subdivisions
    
end

if assembly_input.performance.debug_mode
    abs_error_1 = NaN(matrix_props.n_col, tot_elems);
    rel_error_1 = NaN(matrix_props.n_col, tot_elems);
    fun_evals_1 = NaN(matrix_props.n_col, tot_elems);
    flag_1 =      NaN(matrix_props.n_col, tot_elems);
    min_value_1 = NaN(matrix_props.n_col, tot_elems);
    max_value_1 = NaN(matrix_props.n_col, tot_elems);
    %     time_1 =      NaN(matrix_props.n_col, tot_elems);
    
else
    abs_error_1 = [];
    rel_error_1 = [];
    fun_evals_1 = [];
    flag_1 =      [];
    min_value_1 = [];
    max_value_1 = [];
    %     time_1 = [];
end
if assembly_input.performance.verbose
    disp('initialized Ax Ay Az');   %we can get this far for at least up to 16340 total verts, which should be enough to do dino mesh!  (so far, it has up to 14500 verts)
end
% A = zeros(matrix_props.n_rows, matrix_props.n_cols); %includes extra variables and equations for freeswim case
% disp('initialized A'); %assembling full A will have to be un-mexed for large dino meshes

A_freeswim_rows = zeros(matrix_props.n_rows - matrix_props.n_col*3, matrix_props.n_cols);  %last few rows, all cols
A_freeswim_cols = zeros(matrix_props.n_rows, matrix_props.n_cols - matrix_props.n_col*3); % all rows, last few cols

coder.varsize('VRTS',[2 3 Inf],[true true true]); %for mex rules
VRTS = [0 1 0; 0 0 1];  %reference triangle

rule = 2;  %1st rule has ridiculously innaccurate (too conservative) error estimates, so never converges to the tolerances
[rule_constants.G, rule_constants.Weights, rule_constants.PTS] =  SMPRMS( 2, rule );  %always 2 dimensions and constant rule # for all integrals (2nd arg is rule #)

[vert_range, elem_range] = bounds_from_global_inds(Mesh);  %make it easier to find current submesh inside parfor
[i_mesh_elems, local_elems] = global2local(1:tot_elems, Mesh, elem_range, 'elem');    %compute i_mesh_elem and local_elem, indices to current element's submesh and local element index

if ~assembly_input.skip_traction_integrals %if this is on, we skip the entire core part of A since we can copy it from the previous run (if we are doing freeswim now and just did forced)
    
    for i_mesh_vert = 1:length(Mesh) %only parfor over one submesh at a time to avoid load balancing problem when skip_rigid_integrals = true
        
        col_inds =  Mesh(i_mesh_vert).indices.glob.unq_bounds.vert(1) : Mesh(i_mesh_vert).indices.glob.unq_bounds.vert(2);  %global indices of collocation points for current submesh
        [~, temp] = global2local(col_inds, Mesh, vert_range, 'vert');
        local_verts = NaN(col_inds(end), 1);
        local_verts(col_inds) = temp;  %slot in actual values into larger local_verts vector so we can directly index into it using col_i and avoid parfor errors
        
        mesh_loop_tic = tic;
        
        parfor (col_i = Mesh(i_mesh_vert).indices.glob.unq_bounds.vert(1) : Mesh(i_mesh_vert).indices.glob.unq_bounds.vert(2) , assembly_input.performance.nthreads)  %rows of A, i.e. vertices
            
            %              for col_i = Mesh(i_mesh_vert).indices.glob.unq_bounds.vert(1) : Mesh(i_mesh_vert).indices.glob.unq_bounds.vert(2)
            % col_i
            % for col_i = 1:matrix_props.n_col
            % col_i / matrix_props.n_col
            
            % innertic = tic;
            
            AE = [];   NV = [];     FL = logical([]); %stop parfor complaints about temp variables
           
            
            %             [~, local_vert] = global2local(col_i, Mesh, vert_range, 'vert');
            local_vert = local_verts(col_i);
            
            % [~, local_vert] = global2local(col_i, Mesh, vert_range);     %get local vertex index
            %     ind = matrix_props.global2local.vert.permuted_inds(col_i);  %go from ordered col_i to a randomly permuted global index
            %     i_mesh_vert = matrix_props.global2local.vert.submesh(ind); %random submesh
            %     local_vert = matrix_props.global2local.vert.local(ind);  %random vertex
            
            
            integrand_constants = struct('eps2',assembly_input.accuracy.eps2,'x_col',Mesh(i_mesh_vert).verts(local_vert,:)');
            
            %temporary vectors to later insert into submatrices
            Axtemp = zeros(matrix_props.n_col*3,1); Aytemp = zeros(matrix_props.n_col*3,1); Aztemp = zeros(matrix_props.n_col*3,1);
            
            
            %for debugging + mex rules
            abs_error_temp = [];
            rel_error_temp = [];
            fun_evals_temp = [];
            flag_temp = [];
            min_value_temp = [];
            max_value_temp = [];
            %             time_temp = [];
            
            if assembly_input.performance.debug_mode
                %                 integral = NaN(3,18,tot_elems);  abserror = integral;
                %                 numevals = NaN(1,tot_elems);
                %                 exitflag = numevals;
                
                abs_error_temp = NaN(1,tot_elems);
                rel_error_temp = NaN(1,tot_elems);
                fun_evals_temp = NaN(1,tot_elems);
                flag_temp =      NaN(1,tot_elems);
                min_value_temp = NaN(1,tot_elems);
                max_value_temp = NaN(1,tot_elems);
                %                 time_temp = NaN(1,tot_elems);
                
            end
            
            %             elem_t1 = [];  %appease parfor rules
            
            for elem_i = 1:tot_elems  %elem_i is a global index
                %         elem_i / Mesh.n_elem
                %                 if assembly_input.performance.debug_mode
                %                     elem_t1 = clock;
                %                 end
                
                %                 [i_mesh_elem, local_elem] = global2local(elem_i, Mesh, elem_range, 'elem');    %compute i_mesh_elem and local_elem, indices to current element's submesh and local element index
                i_mesh_elem = i_mesh_elems(elem_i);   local_elem = local_elems(elem_i);
                
                %          ind2 = matrix_props.global2local.elem.permuted_inds(elem_i);  % probably a dummy operation since elem_i = ind2 unless you decide to randomly permute elems as well as verts
                %     i_mesh_elem = matrix_props.global2local.elem.submesh(ind2); %submesh index
                %     local_elem = matrix_props.global2local.elem.local(ind2);  %element index
                
                
                if assembly_input.ignore_interaction && i_mesh_vert ~= i_mesh_elem %ignore interaction between any submeshes
                    
                    %collocation point and element are in different submeshes - skip integrals
                    continue
                end
                
                
                
                
                if assembly_input.skip_rigid_integrals &&  i_mesh_vert == i_mesh_elem %integrals over elements on same object as current collocation point,
                    %which can instead be quickly calculated by rotating values from first timestep
                    %                     continue
                    %collocation point and element are in different submeshes - skip integrals
                    %                             if  assembly_input.performance.debug_mode
                    %                                 disp('skipping integration');
                    %                             end
                    continue
                    
                    %                         elseif  assembly_input.performance.debug_mode
                    %                             disp('not skipping integration');
                    %                         end
                    
                end
                
                %yes, that's right - if both ignore_interaction and
                %skip_rigid_integrals are true, everything is skipped and
                %after the initial rigid integral computation, the code
                %will run fast indeed!
                
                %want i_mesh_elem below since we are generally interested in current
                %element; vertex only matters as far as it's distance away when computing reg stokeslet
                %function
                subverts = Mesh(i_mesh_elem).verts(Mesh(i_mesh_elem).elems(local_elem,:),:);
                
                shape_parameters = Mesh(i_mesh_elem).elem_params(:,local_elem);  %[alpha beta gamma]
                
                col_inds = matrix_props.Col_inds(elem_i,:); %Col_inds was defined based on global indices so use elem_i not local_elem
                
                %scale abs integration tolerance by element area, since big elements
                %should be allowed more total error than small elements
                abstol =  assembly_input.accuracy.integration_tol.traction.abstol * Mesh(i_mesh_elem).area(local_elem);
                %don't scale reltol, since reltol should already account for area
                %since the integral itself should scale with area
                
                %the important part:  adaptively compute boundary integrals of reg stokeslet function over
                %current element, with respect to current collocation point
                
                if assembly_input.performance.debug_mode
                    % shatlab hangs at almost no CPU usage with the tic/toc
                    % inside adsimp so timing inside the parfor doesn't work even though it
                    % compiles without errors
                    %[VL, AE, NV, FL, Time] = adsimp( 2, VRTS, 54,  assembly_input.accuracy.integration_tol.traction.maxevals,abstol, assembly_input.accuracy.integration_tol.traction.reltol ,rule  , rule_constants ,subverts,shape_parameters,integrand_constants,'reg_stokeslet', true);
                    [VL, AE, NV, FL] = adsimp( 2, VRTS, 54,  assembly_input.accuracy.integration_tol.traction.maxevals,abstol, assembly_input.accuracy.integration_tol.traction.reltol ,rule  , rule_constants ,subverts,shape_parameters,integrand_constants,'reg_stokeslet');
                    %                     VL = ones(54,1);  AE = ones(54,1);  NV = 1;  FL = 1;
                else %don't bother storing extra outputs
                    %                     [VL, ~, ~, ~, ~] = adsimp( 2, VRTS, 54,  assembly_input.accuracy.integration_tol.traction.maxevals,abstol, assembly_input.accuracy.integration_tol.traction.reltol ,rule  , rule_constants ,subverts,shape_parameters,integrand_constants,'reg_stokeslet', false);
                    [VL, ~, ~, ~] = adsimp( 2, VRTS, 54,  assembly_input.accuracy.integration_tol.traction.maxevals,abstol, assembly_input.accuracy.integration_tol.traction.reltol ,rule  , rule_constants ,subverts,shape_parameters,integrand_constants,'reg_stokeslet');
                    %                 VL = ones(54,1);  AE = ones(54,1);  NV = 1;  FL = 1;
                    
                end
                
                submat = reshape(VL,18,3);
                
                %entries in matrix are total integrals so can do a cumulative sum
                Axtemp(col_inds) = Axtemp(col_inds) + submat(:,1);
                Aytemp(col_inds) = Aytemp(col_inds) + submat(:,2);
                Aztemp(col_inds) = Aztemp(col_inds) + submat(:,3);
                
                if assembly_input.performance.debug_mode
                    
                    abs_error_temp(elem_i) = max(abs(AE));  %the largest abs error over all 54 component integrals
                    rel_error_temp(elem_i) = max(    abs(AE ./ VL)  );  %the largest rel error over all 54 component integrals
                    fun_evals_temp(elem_i) = NV;
                    flag_temp(elem_i) = FL;
                    
                    min_value_temp(elem_i) = min(abs(VL) );  %smallest (absolute) component integral
                    max_value_temp(elem_i) = max(abs(VL) );  %largest (absolute) component integral
                    
                    %                     time_temp(elem_i) = Time;
                    
                    
                end
                
                
            end  % elements loop
            
            
            numels = matrix_props.n_col * matrix_props.n_col*3;
            %%
            if numels < assembly_input.performance.numels_max  %good enough to just use Ax, Ay, Az
                
                Ax_1(col_i,:) = Axtemp;
                Ay_1(col_i,:) = Aytemp;
                Az_1(col_i,:) = Aztemp;
                
            elseif ceil(numels / 2)  < assembly_input.performance.numels_max  %what if we split Ax, Ay, Az into two, are they small enough then?
                n = ceil(matrix_props.n_col*3 / 2);
                Ax_1(col_i,:) = Axtemp(1:n);
                Ay_1(col_i,:) = Aytemp(1:n);
                Az_1(col_i,:) = Aztemp(1:n);
                
                Ax_2(col_i,:) = Axtemp(n+1 : end);
                Ay_2(col_i,:) = Aytemp(n+1 : end);
                Az_2(col_i,:) = Aztemp(n+1 : end);
                
            else %ceil(matrix_props.n_col * matrix_props.n_col*3 / 4)  < assembly_input.performance.numels_max  %what if we split Ax, Ay, Az into two, are they small enough then?
                n = ceil(matrix_props.n_col*3 / 4);
                
                Ax_1(col_i,:) = Axtemp(1:n);
                Ay_1(col_i,:) = Aytemp(1:n);
                Az_1(col_i,:) = Aztemp(1:n);
                
                Ax_2(col_i,:) = Axtemp(n+1 : n * 2);
                Ay_2(col_i,:) = Aytemp(n+1 : n * 2);
                Az_2(col_i,:) = Aztemp(n+1  :  n * 2);
                
                Ax_3(col_i,:) = Axtemp(n*2+1 : n * 3);
                Ay_3(col_i,:) = Aytemp(n*2+1 : n * 3);
                Az_3(col_i,:) = Aztemp(n*2+1 : n * 3);
                
                Ax_4(col_i,:) = Axtemp(n*3+1 : end);
                Ay_4(col_i,:) = Aytemp(n*3+1 : end);
                Az_4(col_i,:) = Aztemp(n*3+1 : end);
                
            end
            
            %%
            
            if assembly_input.performance.debug_mode
                
                numels =  matrix_props.n_col * tot_elems;
                
                if numels < assembly_input.performance.numels_max
                    
                    abs_error_1(col_i, :) = abs_error_temp(1,:);  %the largest abs error over all 54 component integrals
                    rel_error_1(col_i, :) = rel_error_temp(1,:);  %the largest rel error over all 54 component integrals
                    fun_evals_1(col_i, :) = fun_evals_temp(1,:);
                    flag_1(col_i, :) = flag_temp(1,:);
                    
                    min_value_1(col_i, :) = min_value_temp(1,:);  %smallest (absolute) component integral
                    max_value_1(col_i, :) = max_value_temp(1,:);  %largest (absolute) component integral
                    
                    %                    time_1(col_i,:) = time_temp(1,:);  %time spent for each elem iter
                    
                    % outputting all these bastards will be truly painful,
                    % so leaving unfinished for now unless it really
                    % becomes necessary
                elseif ceil(numels / 2)  < assembly_input.performance.numels_max
                    
                    %                     n = ceil(tot_elems / 2);
                    %                     abs_error_1(col_i, :) = abs_error_temp(1:n);  %the largest abs error over all 54 component integrals
                    %                     rel_error_1(col_i, :) = rel_error_temp(1:n);  %the largest rel error over all 54 component integrals
                    %                     fun_evals_1(col_i, :) = fun_evals_temp(1:n);
                    %                     flag_1(col_i, :) = flag_temp(1:n);
                    %
                    %                     min_value_1(col_i, :) = min_value_temp(1:n);  %smallest (absolute) component integral
                    %                     max_value_1(col_i, :) = max_value_temp(1:n);  %largest (absolute) component integral
                    %
                    %                     abs_error_2(col_i, :) = abs_error_temp(n+1:end);  %the largest abs error over all 54 component integrals
                    %                     rel_error_2(col_i, :) = rel_error_temp(n+1:end);  %the largest rel error over all 54 component integrals
                    %                     fun_evals_2(col_i, :) = fun_evals_temp(n+1:end);
                    %                     flag_2(col_i, :) = flag_temp(n+1:end);
                    %
                    %                     min_value_2(col_i, :) = min_value_temp(n+1:end);  %smallest (absolute) component integral
                    %                     max_value_2(col_i, :) = max_value_temp(n+1:end);  %largest (absolute) component integral
                    
                else
                    %                     n = ceil(tot_elems / 4);
                    %
                    %                     abs_error_1(col_i, :) = abs_error_temp(1:n);  %the largest abs error over all 54 component integrals
                    %                     rel_error_1(col_i, :) = rel_error_temp(1:n);  %the largest rel error over all 54 component integrals
                    %                     fun_evals_1(col_i, :) = fun_evals_temp(1:n);
                    %                     flag_1(col_i, :) = flag_temp(1:n);
                    %
                    %                     min_value_1(col_i, :) = min_value_temp(1:n);  %smallest (absolute) component integral
                    %                     max_value_1(col_i, :) = max_value_temp(1:n);  %largest (absolute) component integral
                    %
                    %                     abs_error_2(col_i, :) = abs_error_temp(n+1:n*2);  %the largest abs error over all 54 component integrals
                    %                     rel_error_2(col_i, :) = rel_error_temp(n+1:n*2);  %the largest rel error over all 54 component integrals
                    %                     fun_evals_2(col_i, :) = fun_evals_temp(n+1:n*2);
                    %                     flag_2(col_i, :) = flag_temp(n+1:n*2);
                    %
                    %                     min_value_2(col_i, :) = min_value_temp(n+1:n*2);  %smallest (absolute) component integral
                    %                     max_value_2(col_i, :) = max_value_temp(n+1:n*2);  %largest (absolute) component integral
                    %
                    %                     abs_error_3(col_i, :) = abs_error_temp(n*2+1 : n * 3);  %the largest abs error over all 54 component integrals
                    %                     rel_error_3(col_i, :) = rel_error_temp(n*2+1 : n * 3);  %the largest rel error over all 54 component integrals
                    %                     fun_evals_3(col_i, :) = fun_evals_temp(n*2+1 : n * 3);
                    %                     flag_3(col_i, :) = flag_temp(n*2+1 : n * 3);
                    %
                    %                     min_value_3(col_i, :) = min_value_temp(n*2+1 : n * 3);  %smallest (absolute) component integral
                    %                     max_value_3(col_i, :) = max_value_temp(n*2+1 : n * 3);  %largest (absolute) component integral
                    %
                    %                     abs_error_4(col_i, :) = abs_error_temp(n*3+1 : end);  %the largest abs error over all 54 component integrals
                    %                     rel_error_4(col_i, :) = rel_error_temp(n*3+1 : end);  %the largest rel error over all 54 component integrals
                    %                     fun_evals_4(col_i, :) = fun_evals_temp(n*3+1 : end);
                    %                     flag_4(col_i, :) = flag_temp(n*3+1 : end);
                    %
                    %                     min_value_4(col_i, :) = min_value_temp(n*3+1 : end);  %smallest (absolute) component integral
                    %                     max_value_4(col_i, :) = max_value_temp(n*3+1 : end);  %largest (absolute) component integral
                    
                end
                
            end
            
            % can display timing debug info, but only works for non-mexed version
            %         task = getCurrentTask;
            %         ID = task.ID;
            %         disp(['Worker = ',num2str(ID),'      ','col_i = ',num2str(col_i),'   took ',num2str(toc(innertic))]);
            
            
        end %parfor over vertices i.e. collocation points %%%%%%%%%%%%%%%%%%%%%%%%%    parfor      %%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        if assembly_input.performance.verbose
            disp(['Matrix assembly for ',Mesh(i_mesh_vert).name,' collocation pts took ',num2str(toc(mesh_loop_tic))]);
        end
    end %for over submeshes
    if assembly_input.performance.verbose
        disp(['Entire matrix assembly loop took ',num2str(toc)]);
    end
    tic;
    
    if assembly_input.performance.debug_mode
        %             Relerror = abs(Abserror ./ Integral);
        %             relcontrol = sum(Relerror(:) <= assembly_input.accuracy.integration_tol.traction.reltol & Abserror(:) >= assembly_input.accuracy.integration_tol.traction.abstol) ./ numel(Abserror);
        %             abscontrol = sum(Relerror(:) >= assembly_input.accuracy.integration_tol.traction.reltol & Abserror(:) <= assembly_input.accuracy.integration_tol.traction.abstol) ./ numel(Abserror);
        %             bothcontrol = sum(Relerror(:) <= assembly_input.accuracy.integration_tol.traction.reltol & Abserror(:) <= assembly_input.accuracy.integration_tol.traction.abstol) ./ numel(Abserror);
        
        
        
        debug_info.abs_error = abs_error_1;  abs_error_1 = [];
        debug_info.rel_error = rel_error_1;  rel_error_1 = [];
        debug_info.fun_evals = fun_evals_1;  fun_evals_1 = [];
        debug_info.flag      = flag_1;       flag_1 = [];
        debug_info.min_value = min_value_1;  min_value_1 = [];
        debug_info.max_value = max_value_1;  max_value_1 = [];
        %         debug_info.time = time_1;            time_1 = [];
        
    end
    
    
    
    %     A(1:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Ax;
    %     %Ax = [];
    %     A(2:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Ay;
    %     %Ay = [];
    %     A(3:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Az;
    %     %Az = [];
    
    
    
    % A = assembly_input.constants.multfactor*A;  %multiply by 1/8/pi as per eq. 3.4, DJ Smith "A boundary element regularized..." paper
    % disp(['Debug info and Ax, Ay, Az building took ',num2str(toc)]);
    
end  %if not skipping all traction integrals


%either way, compute A_force and A_torque, either for direct output and
%later use in forced case, or immediate insertion into A matrix for
%free-swimming case
%tic

[~, A_force]  = compute_force_integral(Mesh, matrix_props, assembly_input); %parallelized   A_force is separated via 3rd dimension for each submesh component
[~, A_torque] = compute_torque_integral(Mesh, matrix_props, assembly_input); %parallelized  A_torque is separated via 3rd dimension for each submesh component
%disp(['A_force and A_torque calculation took ',num2str(toc)])
% switch assembly_input.bugtype
%     case 'dino'
% A_force = rand(3, matrix_props.n_cols - 6,1);  A_torque = A_force;
%     case 'bacteria'
%       A_force = rand(3, matrix_props.n_cols - 6,);  A_torque = A_force;
% end
%A_force = rand(3, matrix_props.n_cols - 6,length(Mesh));  A_torque = A_force;
%%  additional free swimming equations

A_motor_torque = [];  %stop mex complaining
if strcmp(assembly_input.problemtype,'freeswim')
    tic
    
    %     A(matrix_props.n_col*3+1:matrix_props.n_col*3+3,  1:matrix_props.n_col*3) = sum(A_force,3); %sum over submeshes for integrals over entire mesh
    %     A(matrix_props.n_col*3+4:matrix_props.n_col*3+6,  1:matrix_props.n_col*3) = sum(A_torque,3);
    
    A_freeswim_rows(1:3,  1:matrix_props.n_col*3) = sum(A_force,3); %sum over submeshes for integrals over entire mesh
    A_freeswim_rows(4:6,  1:matrix_props.n_col*3) = sum(A_torque,3);
    
    % sum of torques on tail  = specified motor torque (only need for
    % last equation of torque specified BC
    
    %do this no matter what and output row of A below for
    %rotationrate BC, but don't actually add row to A
    
    switch assembly_input.bugtype
        case 'bacteria'
            % get component of tail torque in direction of rotation axis (only for
            % last equation of torque-specified BC)
            tailtorque = A_torque(:,:,2);  %tail hardcoded as 2nd submesh for now but will have to generalize if modeling bacteria with different configuations
            %tailtorque = rand(3,matrix_props.n_cols - 6);
            
            
            switch assembly_input.tail.motorBC
                case 'torque'
                    A_freeswim_rows(7,  1:matrix_props.n_col*3) = motor_orientation(1)*tailtorque(1,:) + motor_orientation(2)*tailtorque(2,:) + motor_orientation(3)*tailtorque(3,:);
                case 'freq'
                    A_motor_torque =  motor_orientation(1)*tailtorque(1,:) + motor_orientation(2)*tailtorque(2,:) + motor_orientation(3)*tailtorque(3,:);
                    %multiplied by f = traction vector, this gives motor torque
            end
    end
    
    temp1 = zeros(matrix_props.n_col, 3); temp2 = temp1;  temp3 = temp1;
    if strcmp(assembly_input.bugtype,'bacteria') && strcmp(assembly_input.tail.motorBC,'torque')
        temp4 = zeros(matrix_props.n_col, 1);  %loop over all n_col for speed, even though verts not on tail will be skipped
        temp5 = temp4;
        temp6 = temp5;
        
    else
        temp4 = zeros(0,1);  temp5 = temp4;  temp6 = temp4;  %stop mex from complaining
    end
    
    [i_mesh_verts, local_verts] = global2local(1:matrix_props.n_col, Mesh, vert_range,'vert');
    
    % move unknown U1 U2 U3 W1 W2 W3 from RHS into A
    parfor (col_i = 1:matrix_props.n_col , assembly_input.performance.nthreads)
        %  for col_i = 1:matrix_props.n_col
        
        i_mesh_vert = i_mesh_verts(col_i);
        local_vert = local_verts(col_i);
        
        %         ind = matrix_props.global2local.vert.permuted_inds(col_i);  %go from ordered col_i to a randomly permuted global index
        %         i_mesh_vert = matrix_props.global2local.vert.submesh(ind); %random submesh
        %         local_vert = matrix_props.global2local.vert.local(ind);  %random vertex
        
        
        x_col = Mesh(i_mesh_vert).verts(local_vert,:)';
        r = x_col - refpoint;
        
        % after n_col*3, columns in A are:  U1 U2 U3 W1 W2 W3 Omega
        %following coeffs for U and W apply for both body and tail
        
        temp1(col_i, :) = -[1 +r(3) -r(2)];
        temp2(col_i, :) = -[1 +r(1) -r(3)];
        temp3(col_i, :) = -[1 +r(2) -r(1)];
        
        if strcmp(assembly_input.bugtype,'bacteria') && strcmp(assembly_input.tail.motorBC,'torque') && i_mesh_vert == 2 %need to solve for omega and col_i is on tail (hardcoded as submesh #2 now, need to generalize for more other bacteria models)
           % - (motor_orientation X r)
            temp4(col_i) = -[motor_orientation(2)*r(3) - motor_orientation(3)*r(2)]; %non-tail verts remain zero and are copied back into A below needlessly, but this method allows the use of just one parfor loop
            temp5(col_i) = -[motor_orientation(3)*r(1) - motor_orientation(1)*r(3)];
            temp6(col_i) = -[motor_orientation(1)*r(2) - motor_orientation(2)*r(1)];
        end
        
    end  %n_col parfor
    
    A_freeswim_cols(3*((1:matrix_props.n_col)-1)+1,  [1 5 6]) = temp1;
    A_freeswim_cols(3*((1:matrix_props.n_col)-1)+2,  [2 6 4]) = temp2;
    A_freeswim_cols(3*((1:matrix_props.n_col)-1)+3, [3 4 5]) = temp3;
    
    if strcmp(assembly_input.bugtype,'bacteria') && strcmp(assembly_input.tail.motorBC,'torque')  %hard coded to assume that Mesh(2) is tail
        A_freeswim_cols(3*( (1:matrix_props.n_col)-1 )+1, 7) = temp4;
        A_freeswim_cols(3*( (1:matrix_props.n_col)-1 )+2,  7) = temp5;
        A_freeswim_cols(3*( (1:matrix_props.n_col)-1 )+3,  7) = temp6;
    end
    
    
    % disp(['Freeswimming part of A took ',num2str(toc)]);
    
    
end  %freeswim


