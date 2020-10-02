function [D] = diffusivities(solutions)

%based on forced translation and rotation flow cases, calculates V, the
%eigenvectors i.e. principal axes of rotation and D, the diagonalized
%rotational diffusivity matrix

kB = 1.3806488E-23 * (1E6 / 1)^2;  %Boltzmann's constant, kg micron^2 / s^2 / K
T = 298;  %temperature to use for diffusivity calc (K)


%Diffusion coefficients for rigid macromolecules with irregular shapes that
%allow rotational-translational coupling


drag.x = [solutions.x.forces.drag];
drag.y = [solutions.y.forces.drag];
drag.z = [solutions.z.forces.drag];

torque.x = [solutions.x.forces.torque];
torque.y = [solutions.y.forces.torque];
torque.z = [solutions.z.forces.torque];

drag.rx = [solutions.rx.forces.drag];
drag.ry = [solutions.ry.forces.drag];
drag.rz = [solutions.rz.forces.drag];

torque.rx = [solutions.rx.forces.torque];
torque.ry = [solutions.ry.forces.torque];
torque.rz = [solutions.rz.forces.torque];



%Kt = translational friction tensor (should be symmetric)
%Kr = rotational friction tensor (should be symmetric)
%Kc = coupling friction tensor
Kt(1,1) = -sum(drag.x(1,:)) / solutions.x.U0;
Kt(2,1) = -sum(drag.x(2,:)) / solutions.x.U0;
Kt(3,1) = -sum(drag.x(3,:)) / solutions.x.U0;
Kc(1,1) = -sum(torque.x(1,:)) / solutions.x.U0;
Kc(2,1) = -sum(torque.x(2,:)) / solutions.x.U0;
Kc(3,1) = -sum(torque.x(3,:)) / solutions.x.U0;

Kt(1,2) = -sum(drag.y(1,:)) / solutions.y.U0;
Kt(2,2) = -sum(drag.y(2,:)) / solutions.y.U0;
Kt(3,2) = -sum(drag.y(3,:)) / solutions.y.U0;
Kc(1,2) = -sum(torque.y(1,:)) / solutions.y.U0;
Kc(2,2) = -sum(torque.y(2,:)) / solutions.y.U0;
Kc(3,2) = -sum(torque.y(3,:)) / solutions.y.U0;

Kt(1,3) = -sum(drag.z(1,:)) / solutions.z.U0;
Kt(2,3) = -sum(drag.z(2,:)) / solutions.z.U0;
Kt(3,3) = -sum(drag.z(3,:)) / solutions.z.U0;
Kc(1,3) = -sum(torque.z(1,:)) / solutions.z.U0;
Kc(2,3) = -sum(torque.z(2,:)) / solutions.z.U0;
Kc(3,3) = -sum(torque.z(3,:)) / solutions.z.U0;

Kr(1,1) = -sum(torque.rx(1,:)) / solutions.rx.Omega0;
Kr(2,1) = -sum(torque.rx(2,:)) / solutions.rx.Omega0;
Kr(3,1) = -sum(torque.rx(3,:)) / solutions.rx.Omega0;

Kr(1,2) = -sum(torque.ry(1,:)) / solutions.ry.Omega0;
Kr(2,2) = -sum(torque.ry(2,:)) / solutions.ry.Omega0;
Kr(3,2) = -sum(torque.ry(3,:)) / solutions.ry.Omega0;

Kr(1,3) = -sum(torque.rz(1,:)) / solutions.rz.Omega0;
Kr(2,3) = -sum(torque.rz(2,:)) / solutions.rz.Omega0;
Kr(3,3) = -sum(torque.rz(3,:)) / solutions.rz.Omega0;

% Kc_check(1,1) = -sum(drag.rx(1,:)) / solutions.rx.Omega0;
% Kc_check(1,2) = -sum(drag.rx(3,:)) / solutions.rx.Omega0;
% Kc_check(1,3) = -sum(drag.rx(3,:)) / solutions.rx.Omega0;
% 
% Kc_check(2,1) = -sum(drag.ry(1,:)) / solutions.ry.Omega0;
% Kc_check(2,2) = -sum(drag.ry(2,:)) / solutions.ry.Omega0;
% Kc_check(2,3) = -sum(drag.ry(3,:)) / solutions.ry.Omega0;
% 
% Kc_check(3,1) = -sum(drag.rz(1,:)) / solutions.rz.Omega0;
% Kc_check(3,2) = -sum(drag.rz(2,:)) / solutions.rz.Omega0;
% Kc_check(3,3) = -sum(drag.rz(3,:)) / solutions.rz.Omega0;


Dt = kB*T*inv(Kt - Kc' * inv(Kr) * Kc);  %translational diffusion matrix (should be symmetric)  eq. 13a
Dr = kB*T*inv(Kr - Kc * inv(Kt) * Kc'); %rotational diffusion matrix (should be symmetric)  eq. 13b
Dc = -inv(Kr) * Kc * Dt; %coupling diffusion matrix  eq. 13c


% %% center of reaction (turns out I don't really need this)
% % eq. 10
% h = [0 0 0]';
% 
% for i = 1:3
%     for j = 1:3
%         for k = 1:3
%             h(i) = h(i) + krondelta(i,j,k) * Kc(j,k);
%         end
%     end
% end
% 
% h = -h;
% 
% cr = inv( trace(Kt)*eye(3) - Kt) * h;  %eq. 9
% %cr = r_OR = center of reaction (for forces) = vector from point P (used to
% %calculate all the friction coeffs) to point R, the center of reaction




%% center of diffusion
% eq 5
d = [0 0 0]';

for i = 1:3
    for j = 1:3
        for k = 1:3
            d(i) = d(i) + krondelta(i,j,k) * Dc(j,k);
        end
    end
end

cd = inv( trace(Dr)*eye(3) - Dr) * d; %eq. 4
%center of diffusion = r_OD = vector from point P (where friction coeffs
%were calculated) to point D, the center of diffusion

%% convert from arbitrarily located values to values at center of diffusion

Dt_cd = Dt - cross2( cross2(cd, Dr), cd) + cross2(Dc', cd) - cross2(cd, Dc); %eq 2a, using cross2 definition from Dave Smith
%translational diffusion matrix evaluated at center of diffusion (should be
%symmetric)

Dc_cd = Dc + cross2(Dr, cd); %eq 2b
%coupling diffusion matrix evaluated at center of diffusion (should be
%symmetric)

% rotational diffusion matrix Dr is apparently independent of position so
% don't need to modify Dr
Dr_cd = Dr;

%% save some intermediate shat
D.details.Kt = Kt;  D.details.Kr = Kr;  D.details.Kc = Kc;
D.details.Dt = Dt;  D.details.Dr = Dr;  D.details.Dc = Dc;
%D.details.Dt_cd = Dt_cd;  D.details.Dr_cd = Dr_cd;  D.details.Dc_cd = Dc_cd;

D.details.presymmetrization.Dt_cd = Dt_cd;
D.details.presymmetrization.Dr_cd = Dr_cd;

%% diagonalize to obtain principle diffusivity values along principle diffusion axes
% 0 for pre-symmetrizing
[D.details.presymmetrization.translation.axes, D.details.presymmetrization.translation.diffusivity] = eig(Dt_cd);
[D.details.presymmetrization.rotation.axes, D.details.presymmetrization.rotation.diffusivity] = eig(Dr_cd);

% since Kt and Kr weren't perfectly symmetric, the axes for both Dt and Dr
% won't be perfectly orthogonal, which is a no-no

% so use Dave's idea to "symmetrize" Dt_cd and Dr_cd and then do eig again

Dt_cd_sym = D.details.presymmetrization.translation.axes * D.details.presymmetrization.translation.diffusivity * D.details.presymmetrization.translation.axes';
Dr_cd_sym = D.details.presymmetrization.rotation.axes * D.details.presymmetrization.rotation.diffusivity * D.details.presymmetrization.rotation.axes';

[D.translation.axes, D.translation.diffusivity] = eig(Dt_cd_sym);
[D.rotation.axes, D.rotation.diffusivity] = eig(Dr_cd_sym);

% even now, if you check dot products between axes vectors, they probably
% won't be exactly zero, but perhaps around 1E-14.  This seems to happen
% with inputs to eig that have non-zero off-diagonal values, even if they
% are exactly symmetric, so there's not much more that can be done to make
% it perfect, but symmetrizing makes the axes much more orthogonal than initially.  
% ('nobalance' option in eig doesn't seem to matter.)


% apparently the order of eigenvalues and eigenvectors is not always x, y,
% z, so make sure it is ordered "correctly"....

% need to update this for pole2pole when tail is often not aligned with x
% axis?



temp  = abs(D.rotation.axes) == repmat(max(abs(D.rotation.axes),[],1),3,1);
[i,j] = ind2sub(size(temp),find(temp == 1));
[i_sort,inds] = sort(i);
D.rotation.axes = D.rotation.axes(:,inds);
D.rotation.diffusivity = D.rotation.diffusivity(:,inds);
D.rotation.diffusivity = diag(D.rotation.diffusivity(D.rotation.diffusivity ~= 0));

temp  = abs(D.translation.axes) == repmat(max(abs(D.translation.axes),[],1),3,1);
[i,j] = ind2sub(size(temp),find(temp == 1));
[i_sort,inds] = sort(i);
D.translation.axes = D.translation.axes(:,inds);
D.translation.diffusivity = D.translation.diffusivity(:,inds);
D.translation.diffusivity = diag(D.translation.diffusivity(D.translation.diffusivity ~= 0));



D.details.Dt_cd = Dt_cd_sym;
D.details.Dr_cd = Dr_cd_sym;


D.coupling.diffusivity = Dc_cd;  %should be all basically zero for non-screwlike shapes.  should always be symmetric

D.center = cd;  %vector from mesh refpoint to center of diffusion




