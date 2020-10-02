function [pt,der, u_n,pt_time ] = sheet_parameterized(t, time, geom)


% def Sinusoid(t):
% 	return [long_lambda*t/(2*math.pi),
% 	long_amp1*math.sin(t - long_omega*time),
% 	0]
%
%
%
% def dSinusoid(t):
% 	return [long_lambda/(2*math.pi),
% 	long_amp1*math.cos(t - long_omega*time),
% 	0]



pt = NaN(length(t),3);  vel = pt;  der = pt;  V_u = pt;

for i = 1:length(t)
    
    %centerline pt
    pt(i,:) = [geom.lambda*t(i)/(2*pi), ...
        ( geom.amp*sin(t(i) - geom.omega*time) ), ...
        0];
    
    % t derivative (tangent dir)  (this is not normalized)
    der(i,:) = [geom.lambda/(2*pi), ...
        geom.amp*cos(t(i) - geom.omega*time), ...
        0];
    
    % normal vector to tangent dir (this is not normalized)
    u_n(i,:) = [ geom.amp*cos(t(i) - geom.omega*time), ...
        - geom.lambda/(2*pi), ...
        0 ];
    
    % time derivative of centerline pt
    pt_time(i,:) = [0, ...
        - geom.omega*geom.amp*cos(t(i) - geom.omega*time), ...
        0];
    
    % velocity of vert is vel of center line pt + vel of vector from centerline
    % pt to vert, so we need to get latter starting from u_n
    
    
    % vector from centerline pt to vert
    
    % u_vert = vert - pt_centerline;
    
    % dist = sqrt(sum(u_vert.^2));  % magnitude of u_vert = distance from centerline to vert
    
    
    % u_n = u_n / sqrt(sum(u_n.^2)) * dist;  %normalize, then scale by dist to make vector from centerline to vert
    % above is no good cause we need analytical u_n to get analytical time
    
    
    
    
end
