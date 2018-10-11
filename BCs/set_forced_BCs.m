function [forced_BCs] = set_forced_BCs(U0,Omega0,flowcase)

%simple function to predefine BCs for certain forced movements

%all rotations are assumed to be around the x, y, z axes going through
%the origin

switch flowcase
    case 'x'  %translation in x
        forced_BCs.U = [U0 0 0]';
        forced_BCs.Omega = [0 0 0]';
    case 'y' %translation in y
        forced_BCs.U = [0 U0 0]';
        forced_BCs.Omega = [0 0 0]';
    case 'z' %translation in z
        forced_BCs.U = [0 0 U0]';
        forced_BCs.Omega = [0 0 0]';
    case 'rx'  %rotation around x
        forced_BCs.U = [0 0 0]';
        forced_BCs.Omega = [Omega0 0 0]';
    case 'ry'  %rotation around y
        forced_BCs.U = [0 0 0]';
        forced_BCs.Omega = [0 Omega0 0]';
    case 'rz'  %rotation around z
        forced_BCs.U = [0 0 0]';
        forced_BCs.Omega = [0 0 Omega0]';
    case 'xrx'  %translation along x + rotation around x
        forced_BCs.U = [U0 0 0]';
        forced_BCs.Omega = [Omega0 0 0]';
        
end

