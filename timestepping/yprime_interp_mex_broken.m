function [der] = yprime_interp_mex(t,y_in,interpolant,nthreads,motor_freq)
% for fast timestepping when interpolated kinematic variables can be used
% use this with ode45 so that time derivatives are taken from trigonmetric
% interpolants, and no further matrix assembly/solve is needed

% motor_freq must be input - use placeholder NaN for a torque BC

if length(y_in) == 6 %constant rotation rate condition, add y(7) internally even though we're not solving for it
    y =  [y_in; t * motor_freq];  %assumes initial angle = 0
else
    y = y_in;
end

y(7) = mod(y(7) + 2*pi, 2*pi);  % shift periodic tail phase angle back to 0 - 2*pi

% c = 0;
% fields = fieldnames(interpolant);
% for f = 1:length(fields)  %each vector variable, e.g. U, Omega, omega
%     for i = 1:length(interpolant.(fields{f}))  %components of each variable, e.g U(1:3), Omega(1:3)
%         c = c + 1;
%         der(c,1) = trig_interp_eval(interpolant.(fields{f})(i),(y(7)));
%     end
% end
%%


der = NaN(length(interpolant),1);
parfor (c = 1:length(interpolant), nthreads)
der(c,1) = trig_interp_eval(interpolant(c),(y(7)));
end

%rotate derivatives from body frame in which they were calculated to
%fixed frame
rotmat = rotation_matrix('z',y(6)) * rotation_matrix('y',y(5)) * rotation_matrix('x',y(4));

der(1:3) = (rotmat * der(1:3));  %convert from body frame to fixed frame
der(4:6) = (rotmat * der(4:6)); %convert from body frame to fixed frame

%do nothing to der(7) which is tail rotation rate, since that is
%always kept in body frame