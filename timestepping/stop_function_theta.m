function [value,isterminal,direction] = stop_function_theta(t,theta)


isterminal = 1; %stop ode45 when value changes sign
direction = 0; %stop regardless of convergence from positive or negative direction

value = 2*pi - theta;  %theta starts at zero and increases monotonically, eventually exceeding 2*pi