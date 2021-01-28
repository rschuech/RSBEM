function [torque] = motor_torque(t, torque0, t_stop)

if t < t_stop
    torque = torque0;
else
    torque = 0;
end