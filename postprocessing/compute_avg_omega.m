function [avg_omega] = compute_avg_omega(interpolant) 

%% compute time-averaged omega over a rotation cycle (and thus for entire trajectory, since omega is in the body frame and thus periodic)
   
    
    fun = @(t,theta) trig_interp_eval(interpolant,theta);
    
    stopfun = @(t,theta) stop_function_theta(t,theta);
    
    sol = ode45(fun,[0 100],0,odeset('InitialStep',1E-12,'events',stopfun,'reltol',1E-10,'abstol',1E-12));
    
%     figure(383)
%     plot(sol.x,sol.y,'o-','markerfacecolor','b'); grid on;  drawnow
    
    T = sol.xe;  %temporal period of a tail rotation (we don't know this a priori since the BC was constant torque, with variable omega)
    avg_omega = 2*pi / T; %avg omega is simply 2*pi radians / time it took to rotate that much
    %%
    
    
   