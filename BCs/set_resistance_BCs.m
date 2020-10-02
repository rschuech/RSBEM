function [resistance_BCs] = set_resistance_BCs(U0,Omega0,flowcase)

%simple function to predefine BCs for certain resistance movements

%all rotations are assumed to be around the x, y, z axes going through some refpoint defined in assemble_RHS

switch flowcase
    case 'x'  %translation in x
        resistance_BCs.U = [U0 0 0]';
        resistance_BCs.Omega = [0 0 0]';
    case 'y' %translation in y
        resistance_BCs.U = [0 U0 0]';
        resistance_BCs.Omega = [0 0 0]';
    case 'z' %translation in z
        resistance_BCs.U = [0 0 U0]';
        resistance_BCs.Omega = [0 0 0]';
    case 'rx'  %rotation around x
        resistance_BCs.U = [0 0 0]';
        resistance_BCs.Omega = [Omega0 0 0]';
    case 'ry'  %rotation around y
        resistance_BCs.U = [0 0 0]';
        resistance_BCs.Omega = [0 Omega0 0]';
    case 'rz'  %rotation around z
        resistance_BCs.U = [0 0 0]';
        resistance_BCs.Omega = [0 0 Omega0]';
    case 'xrx'  %translation along x + rotation around x
        resistance_BCs.U = [U0 0 0]';
        resistance_BCs.Omega = [Omega0 0 0]';
        
end

