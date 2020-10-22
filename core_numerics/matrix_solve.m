function [out1, u, kinematics] = matrix_solve(input, matrix_props, A, RHS)
% outputs f, the solution of the matrix equation
% f includes traction only;
% kinematic values U, Omega, omega are optionally output for free swimming case


matsolve = tic;
%disp(['Solving matrix equation for flowcase ',flowcase]);
sol = A\RHS; %contains both traction f as well as possibly additional unknowns such as U, Omega, omega
if input.performance.verbose
disp(['Matrix solve took ',num2str(toc(matsolve))]);
end
% switch input.problemtype
%     case 'forced'
%  disp(['Matrix solve for flowcase ',flowcase,' took ',num2str(toc(matsolve))]);
%     case 'freeswim'
%          disp(['Matrix solve took ',num2str(toc(matsolve))]);
% end
% now account for value of viscosity  mu, which was assumed to be 1 for
% matrix solve


% since now I use appropriate units so that there shouldn't be any scaling woes, mu is accounting for directly within matrix_assembly_mex
% switch input.problemtype
%     case 'resistance'  %only solved for traction, only need to correct traction
%         sol = sol * input.constants.mu;
%     case 'freeswim'
%         switch input.bugtype
%             case 'bacteria'
%                 switch input.tail.motorBC
%                     case 'freq' %forces change, speeds don't
%                         sol(1:matrix_props.n_col*3) = sol(1:matrix_props.n_col*3) * input.constants.mu;  %correct tractions (get bigger with bigger mu)
%                     case 'torque'  %speeds change, forces don't
%                         sol(matrix_props.n_col*3+1:matrix_props.n_col*3+7) =  sol(matrix_props.n_col*3+1:matrix_props.n_col*3+7) / input.constants.mu;  %correct U, Omega, and freq (get smaller with bigger mu)
%                 end
%             case 'dino'
%                 sol(1:matrix_props.n_col*3) = sol(1:matrix_props.n_col*3) * input.constants.mu;  %correct tractions (get bigger with bigger mu)
%         end
% end

f = sol(1:matrix_props.n_collocation*3 , :);  % only output traction part of sol; kinematic values are output separately if they exist

u = sol(matrix_props.n_collocation*3 + 1: matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3, :); % any unknown free slip velocities (in fixed reference frame)

if input.problemtype == "mobility"
    kinematics.U = sol(matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 1 : matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 3, :); 
    % take all columns in case we've combined several RHS (e.g. 3 translations, 3 rotations in resistance problem)
    kinematics.Omega = sol(matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 4 : matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 6, :); 
    if input.rotating_flagellum && input.Tail.motorBC == "torque"
        kinematics.omega = sol(matrix_props.n_collocation*3 + matrix_props.n_unknown_u*3 + 7, :);
    end
%     varargout = {kinematics};
% else
%     varargout = {};
end

if nargout == 1
    out1 = sol;
else
    out1 = f;
end

%disp(['Solution postprocess took ',num2str(toc(posttic))]);  %takes almost no time