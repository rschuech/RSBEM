function [der] = yprime_interp_mex(t,y_in,interpolant_trig,interpolant_spline,kinematics_interp_method,motor_freq)
% for fast timestepping when interpolated kinematic variables can be used
% use this with ode45 so that time derivatives are taken from trigonmetric
% interpolants, and no further matrix assembly/solve is needed

% motor_freq must be input - use placeholder NaN for a torque BC



coder.extrinsic('fnval');


if length(y_in) == 6 %constant rotation rate condition, add y(7) internally even though we're not solving for it
    y =  [y_in; t * motor_freq];  %assumes initial angle = 0
else
    y = y_in;
end

y(7) = mod(y(7) + 2*pi, 2*pi);  % shift periodic tail phase angle back to 0 - 2*pi

switch kinematics_interp_method
    case 'trig'
        % body frame derivatives
        der_body = NaN(length(interpolant_trig),1);
        % parfor (c = 1:length(interpolant), nthreads)
        for c = 1:length(interpolant_trig)  % slightly faster without parfor but still mexed
            der_body(c,1) = trig_interp_eval(interpolant_trig(c),(y(7)));
        end
    case 'spline'
        der_body = NaN(length(interpolant_spline),1);
        % parfor (c = 1:length(interpolant), nthreads)
        for c = 1:length(interpolant_spline)  % slightly faster without parfor but still mexed
            der_body(c,1) = fnval(interpolant_spline(c),(y(7)));
        end
    otherwise
        der_body = NaN(6,1);
end

% % body frame derivatives
% der_body = NaN(length(interpolant),1);
% % parfor (c = 1:length(interpolant), nthreads)
% for c = 1:length(interpolant)  % slightly faster without parfor but still mexed
%     switch kinematics_interp_method
%         case 'trig'
%             der_body(c,1) = trig_interp_eval(interpolant(c),(y(7)));
%         case 'spline'
%             der_body(c,1) = fnval(interpolant(c),(y(7)));
%     end
% end
der = der_body; %initialization, will overwrite most values


A_1 = A_1_matrix(y(4:6)); % see comments in this function for detailed explanation of what the angles y(4:6) represent

der(1:3) = A_1 * der_body(1:3);  % in the Ramia papers

% see Maple sheet timestepping_eqs_omega
% der(4:6) =   [sin(y(6))/sin(y(5))  ,  cos(y(6))/sin(y(5))  ,  0 ; ...
%     cos(y(6))            ,  -sin(y(6))           ,  0 ; ...
%     -cot(y(5))*sin(y(6)) ,  -cot(y(5))*cos(y(6)) ,  1]  * der_body(4:6);  % in Ramia papers but beware typo/inconsistencies there

der(4:6) = [ 0               ,       sin(y(6))/cos(y(5))  ,       cos(y(6))/cos(y(5))  ; ...
             0               ,       cos(y(6))            ,       -sin(y(6))           ; ...
             1               ,       sin(y(6))*tan(y(5))  ,       cos(y(6))*tan(y(5))  ; ]         * der_body(4:6);
    


% if abs(y(5)) <  1E-10 % 2nd angle is approaching 0 and der will blow up
%     der(4:6) = NaN;
% end

if abs(tan(y(5))) >  1E10 % 2nd angle is approaching 0 and der will blow up
    der(4:6) = NaN;  %hopefully enough to signal a problem
end

% if we're solving for phase angle, it's der got copied from der_body
% already