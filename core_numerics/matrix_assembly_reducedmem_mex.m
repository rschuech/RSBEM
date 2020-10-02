function [Ax,Ay,Az,A_freeswim_rows,A_freeswim_cols,A_motor_torque,A_force,A_torque,debug_info] = matrix_assembly_mex(Mesh,matrix_props,assembly_input)

%assembles all entries of A matrix for both forced and force-free cases
%for free-swim and a freq motorBC , A_motor_torque can be multiplied by f to give the motor
%torque
%for all cases but particularly the forced case, A_force and A_torque can
%be multiplied by f to give the total force and torque on the object
%if debug_mode is true, debug_info can be output


motor_orientation = Mesh(1).orientation(:,1);  % "a" vector for current direction of original x orientation.  Using body but should be same as using tail, since they better both always have same x orientation
refpoint = Mesh(1).refpoints(:,1); %using body but could also use tail as long as refpoint is somewhere along the axis of motor rotation

coder.extrinsic('tic');
coder.extrinsic('toc');
coder.extrinsic('num2str');

debug_info.Integral = [];
debug_info.Abserror = [];
debug_info.Numevals = [];
debug_info.Exitflag = [];
debug_info.Relerror = [];
debug_info.relcontrol = [];
debug_info.abscontrol = [];
debug_info.bothcontrol = [];

tot_elems = 0;
for i = 1:length(Mesh)  %mex is too dumb to understand sum([Mesh.n_elem])...
    tot_elems = Mesh(i).n_elem + tot_elems;
end

tic
%#codegen
%%
%preallocate A matrix and temporary submatrices needed for parallel assembly

Ax = zeros(matrix_props.n_col,matrix_props.n_col*3);  Ay = Ax;  Az = Ax;  %only includes tractions for all cases
disp('initialized Ax Ay Az');   %we can get this far for at least up to 16340 total verts, which should be enough to do dino mesh!  (so far, it has up to 14500 verts)

% A = zeros(matrix_props.n_rows, matrix_props.n_cols); %includes extra variables and equations for freeswim case
% disp('initialized A'); %assembling full A will have to be un-mexed for large dino meshes

A_freeswim_rows = zeros(matrix_props.n_rows - matrix_props.n_col*3, matrix_props.n_cols);  %last few rows, all cols
A_freeswim_cols = zeros(matrix_props.n_rows, matrix_props.n_cols - matrix_props.n_col*3); % all rows, last few cols

%for debug_mode + annoying mex rules
Integral = [];  Abserror = [];  Numevals = [];  Exitflag = [];

if assembly_input.performance.debug_mode
    Integral = NaN(3,18,matrix_props.n_col,tot_elems);
    Abserror = Integral;
    Numevals = NaN(matrix_props.n_col,tot_elems);
    Exitflag = Numevals;
    
end

coder.varsize('VRTS',[2 3 Inf],[true true true]); %for mex rules
VRTS = [0 1 0; 0 0 1];  %reference triangle

rule = 2;  %1st rule has ridiculously innaccurate (too conservative) error estimates, so never converges to the tolerances
[rule_constants.G, rule_constants.Weights, rule_constants.PTS] =  SMPRMS( 2, rule );  %always 2 dimensions and constant rule # for all integrals (2nd arg is rule #)

[vert_range, elem_range] = bounds_from_global_inds(Mesh);  %make it easier to find current submesh inside parfor

if ~assembly_input.skip_traction_integrals %if this is on, we skip the entire core part of A since we can copy it from the previous run (if we are doing freeswim now and just did forced)
    
    for i_mesh_vert = 1:length(Mesh) %only parfor over one submesh at a time to avoid load balancing problem when skip_rigid_integrals = true
        
        parfor (col_i = Mesh(i_mesh_vert).global_indices.vert.start : Mesh(i_mesh_vert).global_indices.vert.end , assembly_input.performance.nthreads)  %rows of A, i.e. vertices
            % for col_i = 1:matrix_props.n_col
            % col_i / matrix_props.n_col
            
            % innertic = tic;
            
            AE = [];   NV = []; FL = logical([]); %stop parfor complaints about temp variables
            
            [~, local_vert] = global2local(col_i, Mesh, vert_range);     %get local vertex index
            %     ind = matrix_props.global2local.vert.permuted_inds(col_i);  %go from ordered col_i to a randomly permuted global index
            %     i_mesh_vert = matrix_props.global2local.vert.submesh(ind); %random submesh
            %     local_vert = matrix_props.global2local.vert.local(ind);  %random vertex
            
            
            integrand_constants = struct('eps2',assembly_input.accuracy.eps2,'x_col',Mesh(i_mesh_vert).verts(local_vert,:)');
            
            %temporary vectors to later insert into submatrices
            Axtemp = zeros(1,matrix_props.n_col*3);  Aytemp = Axtemp;  Aztemp = Axtemp;
            
            
            %for debugging + mex rules
            integral = [];  abserror = [];  numevals = [];  exitflag = [];
            
            if assembly_input.performance.debug_mode
                integral = NaN(3,18,tot_elems);  abserror = integral;
                numevals = NaN(1,tot_elems);
                exitflag = numevals;
            end
            
            for elem_i = 1:tot_elems  %elem_i is a global index
                %         elem_i / Mesh.n_elem
                
                [i_mesh_elem, local_elem] = global2local(elem_i, Mesh, elem_range);    %compute i_mesh_elem and local_elem, indices to current element's submesh and local element index
                %          ind2 = matrix_props.global2local.elem.permuted_inds(elem_i);  % probably a dummy operation since elem_i = ind2 unless you decide to randomly permute elems as well as verts
                %     i_mesh_elem = matrix_props.global2local.elem.submesh(ind2); %submesh index
                %     local_elem = matrix_props.global2local.elem.local(ind2);  %element index
                
                if strcmp(assembly_input.problemtype,'freeswim')
                    if assembly_input.ignore_interaction  %ignore interaction between any submeshes
                        
                        if i_mesh_vert ~= i_mesh_elem %collocation point and element are in different submeshes - skip integrals
                            continue
                        end
                        
                    end
                    
                    if assembly_input.skip_rigid_integrals %integrals over elements on same object as current collocation point,
                        %which can instead be quickly calculated by rotating values from first timestep
                        
                        if i_mesh_vert == i_mesh_elem %collocation point and element are in different submeshes - skip integrals
                            continue
                        end
                        
                    end
                    
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
                    [VL, AE, NV, FL] = adsimp( 2, VRTS, 54,  assembly_input.accuracy.integration_tol.traction.maxevals,abstol, assembly_input.accuracy.integration_tol.traction.reltol ,rule  , rule_constants ,subverts,shape_parameters,integrand_constants,'reg_stokeslet');
                else %don't bother storing extra outputs
                    [VL, ~, ~, ~] = adsimp( 2, VRTS, 54,  assembly_input.accuracy.integration_tol.traction.maxevals,abstol, assembly_input.accuracy.integration_tol.traction.reltol ,rule  , rule_constants ,subverts,shape_parameters,integrand_constants,'reg_stokeslet');
                end
                
                submat = reshape(VL,18,3)';
                
                %entries in matrix are total integrals so can do a cumulative sum
                Axtemp(col_inds) = Axtemp(col_inds) + submat(1,:);
                Aytemp(col_inds) = Aytemp(col_inds) + submat(2,:);
                Aztemp(col_inds) = Aztemp(col_inds) + submat(3,:);
                
                if assembly_input.performance.debug_mode
                    integral(:,:,elem_i) = submat;
                    AE = reshape(AE,18,3)';
                    abserror(:,:,elem_i) = AE;
                    numevals(elem_i) = NV;
                    exitflag(elem_i) = FL;
                end
                
            end  % elements loop
            
            Ax(col_i,:) = Axtemp;
            Ay(col_i,:) = Aytemp;
            Az(col_i,:) = Aztemp;
            
            if assembly_input.performance.debug_mode
                Integral(:,:,col_i,:) = integral;
                Abserror(:,:,col_i,:) = abserror;
                Numevals(col_i,:) = numevals(1,:);
                Exitflag(col_i,:) = exitflag(1,:);
            end
            
            % can display timing debug info, but only works for non-mexed version
            %         task = getCurrentTask;
            %         ID = task.ID;
            %         disp(['Worker = ',num2str(ID),'      ','col_i = ',num2str(col_i),'   took ',num2str(toc(innertic))]);
            
            
        end %parfor over vertices i.e. collocation points %%%%%%%%%%%%%%%%%%%%%%%%%    parfor      %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end %for over submeshes
    
    disp(['Parallel matrix assembly loop took ',num2str(toc)]);
    
    tic;
    
    if assembly_input.performance.debug_mode
        Relerror = abs(Abserror ./ Integral);
        relcontrol = sum(Relerror(:) <= assembly_input.accuracy.integration_tol.traction.reltol & Abserror(:) >= assembly_input.accuracy.integration_tol.traction.abstol) ./ numel(Abserror);
        abscontrol = sum(Relerror(:) >= assembly_input.accuracy.integration_tol.traction.reltol & Abserror(:) <= assembly_input.accuracy.integration_tol.traction.abstol) ./ numel(Abserror);
        bothcontrol = sum(Relerror(:) <= assembly_input.accuracy.integration_tol.traction.reltol & Abserror(:) <= assembly_input.accuracy.integration_tol.traction.abstol) ./ numel(Abserror);
        
        
        
        debug_info.Integral = Integral;
        debug_info.Abserror = Abserror;
        debug_info.Numevals = Numevals;
        debug_info.Exitflag = Exitflag;
        debug_info.Relerror = Relerror;
        debug_info.relcontrol = relcontrol;
        debug_info.abscontrol = abscontrol;
        debug_info.bothcontrol = bothcontrol;
    end
    
    
    
%     A(1:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Ax;
%     %Ax = [];
%     A(2:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Ay;
%     %Ay = [];
%     A(3:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Az;
%     %Az = [];
    
    
    
   % A = assembly_input.constants.multfactor*A;  %multiply by 1/8/pi as per eq. 3.4, DJ Smith "A boundary element regularized..." paper
    disp(['Debug info and A building took ',num2str(toc)]);
    
end  %if not skipping all traction integrals


%either way, compute A_force and A_torque, either for direct output and
%later use in forced case, or immediate insertion into A matrix for
%free-swimming case
tic
[~, A_force] = compute_force_integral(Mesh, matrix_props, assembly_input); %parallelized   A_force is separated via 3rd dimension for each submesh component
[~, A_torque] = compute_torque_integral(Mesh, matrix_props, assembly_input); %parallelized  A_torque is separated via 3rd dimension for each submesh component
disp(['A_force and A_torque calculation took ',num2str(toc)]);

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
    
    
    % get component of tail torque in direction of rotation axis (only for
    % last equation of torque-specified BC)
    tailtorque = A_torque(:,:,2);  %tail hardcoded as 2nd submesh for now but will have to generalize if modeling bacteria with different configuations
    
    switch assembly_input.tail.motorBC
        case 'torque'
            A_freeswim_rows(7,  1:matrix_props.n_col*3) = motor_orientation(1)*tailtorque(1,:) + motor_orientation(2)*tailtorque(2,:) + motor_orientation(3)*tailtorque(3,:);
        case 'freq'
            A_motor_torque =  motor_orientation(1)*tailtorque(1,:) + motor_orientation(2)*tailtorque(2,:) + motor_orientation(3)*tailtorque(3,:);
            %multiplied by f = traction vector, this gives motor torque
    end
    
    
    temp1 = zeros(matrix_props.n_col, 3); temp2 = temp1;  temp3 = temp1;
    if strcmp(assembly_input.tail.motorBC,'torque')
        temp4 = zeros(matrix_props.n_col, 1);  %loop over all n_col for speed, even though verts not on tail will be skipped
        temp5 = temp4;
        temp6 = temp5;
        
    else
        temp4 = zeros(0,1);  temp5 = temp4;  temp6 = temp4;  %stop mex from complaining
    end
    
    % move unknown U1 U2 U3 W1 W2 W3 from RHS into A
    parfor (col_i = 1:matrix_props.n_col , assembly_input.performance.nthreads)
        %  for col_i = 1:matrix_props.n_col
        
        [i_mesh_vert, local_vert] = global2local(col_i, Mesh, vert_range);
        
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
        
        if strcmp(assembly_input.tail.motorBC,'torque') && i_mesh_vert == 2 %need to solve for omega and col_i is on tail (hardcoded as submesh #2 now, need to generalize for more other bacteria models)
            temp4(col_i) = -[motor_orientation(2)*r(3) - motor_orientation(3)*r(2)]; %non-tail verts remain zero and are copied back into A below needlessly, but this method allows the use of just one parfor loop
            temp5(col_i) = -[motor_orientation(3)*r(1) - motor_orientation(1)*r(3)];
            temp6(col_i) = -[motor_orientation(1)*r(2) - motor_orientation(2)*r(1)];
        end
        
    end  %n_col parfor
    
    A_freeswim_cols(3*((1:matrix_props.n_col)-1)+1,  [1 5 6]) = temp1;
    A_freeswim_cols(3*((1:matrix_props.n_col)-1)+2,  [2 6 4]) = temp2;
    A_freeswim_cols(3*((1:matrix_props.n_col)-1)+3, [3 4 5]) = temp3;
    
    if strcmp(assembly_input.tail.motorBC,'torque')  %hard coded to assume that Mesh(2) is tail
        A_freeswim_cols(3*( (1:matrix_props.n_col)-1 )+1, 7) = temp4;
        A_freeswim_cols(3*( (1:matrix_props.n_col)-1 )+2,  7) = temp5;
        A_freeswim_cols(3*( (1:matrix_props.n_col)-1 )+3,  7) = temp6;
    end
    
    
    disp(['Freeswimming part of A took ',num2str(toc)]);
    
    
end  %freeswim


