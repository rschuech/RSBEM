function [V,D,cr] = rotational_diffusivities(solutions,kB,T)

%based on forced translation and rotation flow cases, calculates V, the
%eigenvectors i.e. principal axes of rotation and D, the diagonalized
%rotational diffusivity matrix

%Diffusion coefficients for rigid macromolecules with irregular shapes that
%allow rotational-translational coupling

%Kt = translational friction tensor (should be symmetric)
%Kr = rotational friction tensor (should be symmetric)
%Kc = coupling friction tensor
Kt(1,1) = -solutions.x.forces.drag(1) / solutions.x.U0;
Kt(2,1) = -solutions.x.forces.drag(2) / solutions.x.U0;
Kt(3,1) = -solutions.x.forces.drag(3) / solutions.x.U0;
Kc(1,1) = -solutions.x.forces.torque(1) / solutions.x.U0;
Kc(2,1) = -solutions.x.forces.torque(2) / solutions.x.U0;
Kc(3,1) = -solutions.x.forces.torque(3) / solutions.x.U0;

Kt(1,2) = -solutions.y.forces.drag(1) / solutions.y.U0;
Kt(2,2) = -solutions.y.forces.drag(2) / solutions.y.U0;
Kt(3,2) = -solutions.y.forces.drag(3) / solutions.y.U0;
Kc(1,2) = -solutions.y.forces.torque(1) / solutions.y.U0;
Kc(2,2) = -solutions.y.forces.torque(2) / solutions.y.U0;
Kc(3,2) = -solutions.y.forces.torque(3) / solutions.y.U0;

Kt(1,3) = -solutions.z.forces.drag(1) / solutions.z.U0;
Kt(2,3) = -solutions.z.forces.drag(2) / solutions.z.U0;
Kt(3,3) = -solutions.z.forces.drag(3) / solutions.z.U0;
Kc(1,3) = -solutions.z.forces.torque(1) / solutions.z.U0;
Kc(2,3) = -solutions.z.forces.torque(2) / solutions.z.U0;
Kc(3,3) = -solutions.z.forces.torque(3) / solutions.z.U0;

Kr(1,1) = -solutions.rx.forces.torque(1) / solutions.rx.Omega0;
Kr(2,1) = -solutions.rx.forces.torque(2) / solutions.rx.Omega0;
Kr(3,1) = -solutions.rx.forces.torque(3) / solutions.rx.Omega0;

Kr(1,2) = -solutions.ry.forces.torque(1) / solutions.ry.Omega0;
Kr(2,2) = -solutions.ry.forces.torque(2) / solutions.ry.Omega0;
Kr(3,2) = -solutions.ry.forces.torque(3) / solutions.ry.Omega0;

Kr(1,3) = -solutions.rz.forces.torque(1) / solutions.rz.Omega0;
Kr(2,3) = -solutions.rz.forces.torque(2) / solutions.rz.Omega0;
Kr(3,3) = -solutions.rz.forces.torque(3) / solutions.rz.Omega0;

Kc_check(1,1) = -solutions.rx.forces.drag(1) / solutions.rx.Omega0;
Kc_check(1,2) = -solutions.rx.forces.drag(2) / solutions.rx.Omega0;
Kc_check(1,3) = -solutions.rx.forces.drag(3) / solutions.rx.Omega0;

Kc_check(2,1) = -solutions.ry.forces.drag(1) / solutions.ry.Omega0;
Kc_check(2,2) = -solutions.ry.forces.drag(2) / solutions.ry.Omega0;
Kc_check(2,3) = -solutions.ry.forces.drag(3) / solutions.ry.Omega0;

Kc_check(3,1) = -solutions.rz.forces.drag(1) / solutions.rz.Omega0;
Kc_check(3,2) = -solutions.rz.forces.drag(2) / solutions.rz.Omega0;
Kc_check(3,3) = -solutions.rz.forces.drag(3) / solutions.rz.Omega0;


Dt = kB*T*inv(Kt - Kc' * inv(Kr) * Kc);  %translational diffusion tensor (should be symmetric)  eq. 13a
Dr = kB*T*inv(Kr - Kc * inv(Kt) * Kc'); %rotational diffusion tensor (should be symmetric)  eq. 13b
Dc = -inv(Kr) * Kc * Dt; %coupling difufsion tensor  eq. 13c

%center of reaction

%% eq. 10
h = [0 0 0]';

for i = 1:3
    for j = 1:3
        for k = 1:3
            h(i) = h(i) + krondelta(i,j,k) * Kc(j,k);
        end
    end
end

h = -h;
%%

cr = inv( trace(Kt)*eye(3) - Kt) * h;  %eq. 9
%cr = r_OR = center of reaction (for forces) = vector from point P (used to
%calculate all the friction coeffs) to point R, the center of reaction


%% eq 5
d = [0 0 0]';

for i = 1:3
    for j = 1:3
        for k = 1:3
            d(i) = d(i) + krondelta(i,j,k) * Dc(j,k);
        end
    end
end

%%

cd = inv( trace(Dr)*eye(3) - Dr) * d; %eq. 4
%center of diffusion = r_OD = vector from point P (where friction coeffs
%were calculated) to point D, the center of diffusion

Dt_d = Dt - cross2( cross2(cd, Dr), cd) + cross2(Dc', cd) - cross2(cd, Dc); %eq 2a, using cross2 definition from Dave Smith
%translational diffusion matrix evaluated at center of diffusion


[V,D] = eig(Dr);

stop