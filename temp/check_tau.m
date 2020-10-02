
 swimming_axis = [1 0 0]';  swimming_axis = swimming_axis / sqrt(sum(swimming_axis.^2));
%  D.rotation.axes = [1 0 0; 0 1 0; 0 0 1;];
%  D.rotation.diffusivity = diag([2 3 4]);
%% 
swimming_axis = fits.avg_swimming_axis;

 principle_ind = 3;
 principle_axis_1 =  D.rotation.axes(:,principle_ind); 
% 
% 
%  inds = [1 2];
 inds = setdiff(1:3,principle_ind);
% tau_theory = 1/(  D.rotation.diffusivity(inds(1),inds(1)) + D.rotation.diffusivity(inds(2),inds(2))  )
% 
% 
% angle1 = 90*pi/180;
%  rotmat1 = rotate_arbitrary_vector( [1 0 0]', angle1);
%  angle2 = 90*pi/180;
%   rotmat2 = rotate_arbitrary_vector( [0 1 0]', angle2);
% rotmat = rotmat2*rotmat1;
% 
% axes2 = rotmat * D.rotation.axes
% diff2 = rotmat * D.rotation.diffusivity * rotmat'
% 
% [vecs,vals] = eig(diff2)

%


    rotvec = crossprod( principle_axis_1 , swimming_axis);  rotvec = rotvec / sqrt(sum(rotvec.^2));
    angle = acos( dot(swimming_axis, principle_axis_1) );
    rotmat = rotate_arbitrary_vector( rotvec, angle);
    D_aligned = rotmat * D.rotation.diffusivity * rotmat'  % that's how you rotate a tensor
    axes_aligned = rotmat * D.rotation.axes;
    
%     inds0 = [2 3];
%     1/(  D.rotation.diffusivity(inds0(1),inds0(1)) + D.rotation.diffusivity(inds0(2),inds0(2))  )
%      1/(  D_aligned(inds(1),inds(1)) + D_aligned(inds(2),inds(2))  )
    %
    
 fun = @(angle) - rotate_off_axis_tau(angle, swimming_axis, D_aligned , principle_ind);

guesses = linspace(0,2*pi,10);  angle = NaN(size(guesses));  max_tau = angle;  % 10 evenly spaced initial guesses seems quite enough (3 would probably do)
for nn = 1:length(guesses)
[angle(nn),temp] = fminsearch(fun,guesses(nn));
max_tau(nn) = -temp;
end
[max_tau,ind] = max(max_tau);
angle = angle(ind);

  
    max_tau;
   
    clear tau
    angles = linspace(0,2*pi,500);
    for an = 1:length(angles)
    tau(an) = -fun(angles(an));
    end
    
    figure(305)
    plot(angles,tau,'k-','linewidth',2);  grid on;
    