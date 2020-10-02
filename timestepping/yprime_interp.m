function [der] = yprime_interp(t,y_in,interpolant,kinematics_interp_method,periodic_kinematics,motor_freq)
% for fast timestepping when interpolated kinematic variables can be used
% use this with ode45 so that time derivatives are taken from trigonmetric
% interpolants, and no further matrix assembly/solve is needed

% motor_freq must be input - use placeholder NaN for a torque BC

% turns out that for trig interp, only a very slight speedup by mexing.
% and for spline, fnval isn't mexable and that's where all the time is
% spent (spline is indeed much slower than trig, but still insignificant
% compared to rest of code)
% so not much point in bothering with mexed version


if length(y_in) == 6 %constant rotation rate condition, add y(7) internally even though we're not solving for it
    y =  [y_in; t * motor_freq];  %assumes initial angle = 0.  This is a monotonically increasing phase angle, not yet wrapped to 0 - 2*pi
else
    y = y_in;
end

if ~periodic_kinematics % currently this means we're doing a dino turning simulation; input variable interpolant is actually a piecewise_interpolant struct
    
    ind = find( y(7) >= interpolant.phase_breaks(:,1) & y(7) <= interpolant.phase_breaks(:,2) );
    if isempty(ind) || ~isscalar(ind)
        error('problem with piecewise interpolant definition');
    else
%         if interpolant.interpolant_indices(ind) ~= 1
%             stopafra
%         end
        interpolant = interpolant.interpolants( interpolant.interpolant_indices(ind)  , :);
    end
end

y(7) = mod(y(7) + 2*pi, 2*pi);  % shift periodic tail phase angle back to 0 - 2*pi



% body frame derivatives
der_body = NaN(length(interpolant),1);
% parfor (c = 1:length(interpolant), length(interpolant))  % ubelievably slower than serial!
for c = 1:length(interpolant)  
    switch kinematics_interp_method
        case 'trig'
            der_body(c,1) = trig_interp_eval(interpolant(c),(y(7)));
        case 'spline'
            der_body(c,1) = fnval(interpolant(c),(y(7)));
    end
end
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