function [der,interp_partial] = yprime(t,y,y0,Mesh,matrix_props,assembly_input,BCs,A0,A_motor_torque0, A_recycle, varargin)
%y(1:3) is current refpoint
global solutions
%y(4 5 6) = x, y, z rotation angles of body

%y(7) = rotation angle of tail

% body orientation described by vectors a and b, originally along the x and
% y axes at t = 0.

%% translate and rotate Mesh to current position and orientation

if isempty(varargin)
    
    if length(y) == 6 %constant rotation rate condition, calculate current tail phase angle here if it was not input (i.e. full timestepping simulation)
        switch assembly_input.bugtype
            case 'bacteria'
                y(7) =  t * BCs.motor_freq;  %assumes initial tail phase angle = 0
            case {'dino', 'sheet'}
                y(7) = t * BCs.freeswim.phase_speed;
        end
        % this extra value just added to y is not retained between ode45 calls - it's
        % not part of the ODE solution, but is used for
        % convenience inside this function
    end
    
    
    if strcmp(assembly_input.bugtype,'bacteria') && isequal(y,y0)
        
        A = A_recycle;  %just use originally calculated A, it won't change at all if y hasn't
        A_motor_torque = A_motor_torque0;
        debug_info = [];  %placeholder if we are saving debug_info dumps
        
    else
        temp = tic;
        [A, A_motor_torque, ~, ~, debug_info] = matrix_assembly_mex_wrapper(Mesh,matrix_props,assembly_input);  %rigid integrals might be skipped depending on assembly_input
        if assembly_input.performance.verbose
            disp(['matrix_assembly_mex_wrapper took ',num2str(toc(temp))]);
        end
    end
    
    
    if assembly_input.skip_rigid_integrals  && ~isequal(y,y0)   %rigid integrals (i.e., collocation points and elements on same rigid object) were skipped, so need to fill them in based on rotation of A0
        %A0 has to be orig t = 0 orientation that current orientations in y are in reference to
        rotate_A_tic = tic;
        
        %%
        % instead of recomputing all of A at each time / configuration, can simply rotate the entries
        % corresponding to collocation points and elements on the same submesh, as
        % long as the submesh appears to rigidly rotate over time
        
        %leaves values alone for collocation point and element on different submeshes
        %rotates values for collocation point and element on the input submesh to
        %match current orientation angles
        
        
        % this used to be a function, rotate_A, but that leads to extra copies of A
        % or at least the A0 parts of A, so having everything direclty here avoids
        % any unnecessary copies at expense of ugliness
        
        for ii = 1:length(Mesh)  % do body, then tail, tensor rotations
            switch ii
                case 1
                    angles = y(4:6);  %rotate body entries
                case 2
                    angles = [y(4) y(5) y(6)+y(7)];  %rotate tail entries; recall that tail rotation in x is combo of body rotation and tail phase angle
                    % note update here, since we used to have y(4:7) = [X Y Z phase] but now
                    % y(4:7) = [Z Y X phase] to follow Goldstein notation
            end
            
            [rotmat] = A_1_matrix(angles); % using intrinsic z-y-x rotation convention in Goldstein text
            
            %         rotmat = rotation_matrix('z',angles(3)) * rotation_matrix('y',angles(2)) * rotation_matrix('x',angles(1));
            blk = kron(speye(Mesh(ii).n_vert),rotmat);  %using sparse is *a lot* faster than not - I checked.  also, mexed version of this function is slower - turns out shatlab is actually pretty good at some things
            
            A( (Mesh(ii).indices.glob.unq_bounds.vert(1) - 1)*3+1  :  Mesh(ii).indices.glob.unq_bounds.vert(2) * 3 , ...
                (Mesh(ii).indices.glob.unq_bounds.vert(1) - 1)*3+1  :  Mesh(ii).indices.glob.unq_bounds.vert(2) * 3) ...
                = blk *    A0{ii}    * blk';
            
        end
        %%
        if assembly_input.performance.verbose
            disp(['Tensor rotations took ',num2str(toc(rotate_A_tic))]);
        end
    end
    
    
    
    
    %% given current Mesh, solve for U, Omega, and possibly omega
    [RHS] = assemble_RHS(assembly_input,Mesh, matrix_props, BCs);
    [f, kinematics] = matrix_solve(assembly_input, matrix_props, A, RHS); % traction is f vector and kinematics contains U, Omega, omega
    interp_partial.f = f;  interp_partial.kinematics = kinematics; % send out for saving to file
else
    interp_partial = varargin{1};
    f = interp_partial.f;
    kinematics = interp_partial.kinematics;
    
end

%%

% der(1) = dx/dt = U(1)
% der(2) = dy/dt = U(2)
% der(3) = dz/dt = U(3)

% der(4) = d(theta_x)/dt = Omega(1)
% der(5) = d(theta_y)/dt = Omega(2)
% der(6) = d(theta_z)/dt = Omega(3)

% der(7) = d(theta_tail_phase)/dt = omega    for torque BC case

der(1:3,1) = kinematics.U;

der(4:6,1) = kinematics.Omega;

if strcmp(assembly_input.bugtype,'bacteria') && strcmp(assembly_input.tail.motorBC,'torque')
    der(7,1) = kinematics.omega;
end


%% save outputs besides position:  traction, kinematics, and possibly tail torque
chunksize = 100;
% grow solutions in chunks to mitigate slowdown due to growing arrays
% (should only happen during full timestepping simulations, not
% phase interpolations, for which solutions is repeatedly initialized to correct size = length(phase angles))
if solutions.n_out + 1 > size(solutions.f,1) %if index for next entry is greater than current size of solutions
    fields = fieldnames(solutions);
    fields = setdiff(fields,{'n_out', 'f'}); %n_out is just a scalar, don't grow its size
    for ff = 1:length(fields)
        solutions.(fields{ff})(end+1:end+chunksize,:) = NaN;
    end
    solutions.f{end+1:end+chunksize,1} = [];
end

solutions.n_out = solutions.n_out + 1;

solutions.f{solutions.n_out,1} = f;
if strcmp(assembly_input.bugtype,'bacteria') && strcmp(assembly_input.tail.motorBC,'freq') %torque varies, so store it
    solutions.motor_torque(solutions.n_out, 1) = A_motor_torque * f;
end

fields = fieldnames(kinematics);  % U, Omega, maybe omega
for f = 1:length(fields)
    solutions.(fields{f})(solutions.n_out,:) = kinematics.(fields{f});
end


% switch assembly_input.motorBC
%     case 'torque'
%         solf(n_out,:) = [t; f; tailtorque]; %has rotationrate in f
%     case 'rotationrate'
%         solf(n_out,:) = [t; f; motor_freq; tailtorque];
% end
% elapsed = toc(timerval);
% fraction_done = (t - interval(1)) /  (interval(2) - interval(1));
% fraction_left = 1 - fraction_done;
% timeleft = elapsed / fraction_done * fraction_left;
% disp(['n_out = ',num2str(n_out),'       t = ',num2str(t),'       elapsed = ',num2str(elapsed/60,2),' min','       timeleft = ',num2str(timeleft/60,2),' min']);
% n_out = n_out + 1;
%
% drawnow
% if assembly_input.performance.debug_mode
%     debug_dump_filename = varargin{1};
%     save(debug_dump_filename,'debug_info');
% end





