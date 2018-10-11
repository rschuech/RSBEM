function [fcoeffs] = friction_coeffs(forces,BCs)


% sum the forces, torques on each submesh when computing friction coeffs
% for the entire mesh
fcoeffs.translation = - sum([forces.drag],2) ./ BCs.U;  %three translational friction coeffs for x, y, z (note that only one corresponding to current flowcase will be important)
fcoeffs.rotation = - sum([forces.torque],2) ./ BCs.Omega;  %three rotational friction coeffs for x, y, z (note that only one corresponding to current flowcase will be important)
