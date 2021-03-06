temp = load('C:\Users\rudi\Desktop\RD\Results\body only\Results_body_only_fixed.mat');

Results_body = temp.Results;


 body_geoms = [ [Results_body.AR1]'   [Results_body.AR2]'    ];
 
 
for results_ind = 1:length(Results)

%  Results(results_ind).error_angle = acos(dot( Results(results_ind).principle_axis_1, swimming_axis)) * 180/pi;  %how big is erroneous angle between rotational diffusion principle axis and avg swimming axis
    
    geom_current = [   Results(results_ind).AR1   Results(results_ind).AR2    ];
    
    [~,body_ind] =      ismember(geom_current, body_geoms, 'rows');
    if body_ind == 0
        Results(results_ind).tau_a_body = NaN;
        disp(['Need to run body only forced simulation for AR1 ',num2str(geom_current(1)),'  AR2 ',num2str(geom_current(2))]);
        continue
    end
    
    swimming_axis = Results(results_ind).path_slope ./ sqrt(sum(Results(results_ind).path_slope.^2));
    % using actual swimming axis shifts fore-aft S/N optimum slightly up to
    % AR2 = 0.05 but it should really be straight, turns out that assuming
    % swimming along x axis fixes it.  Largest deviation from x axis is a
    % few degrees so not a big deal.
%     swimming_axis = [1 0 0]';
    
    if isequal(geom_current,[1 0]) %if geom is an exact sphere, need to hardcode principle axis since code outputs a random vector due to symmetry!
        principle_axis_1 = [1 0 0]';
    else
        principle_axis_1 = Results_body(body_ind).principle_axis_1;
    end
    rotvec = crossprod( principle_axis_1 , swimming_axis);  rotvec = rotvec / sqrt(sum(rotvec.^2));
    angle = acos( dot(swimming_axis, principle_axis_1) );
    rotmat = rotate_arbitrary_vector( rotvec, angle);
    D_aligned = rotmat * Results_body(body_ind).rotational_diffusivity * rotmat';  % that's how you rotate a tensor
% turns out that doing following averaging causes a sharp rift in tau for
% unknown reasons, must be due to screwlike properties.  

%     fun = @(angle) rotate_off_axis_tau(angle, swimming_axis, D_aligned, Results_body(body_ind).principle_axis_1_index);
%     limits = [0 2*pi];
    
%     tau_mean = 1/diff(limits) * integral(fun, limits(1),limits(2),'reltol',1E-9,'abstol',1E-12);
%       rotmat = rotate_arbitrary_vector(swimming_axis, pi);
%     
%     D_rot = rotmat * D_aligned * rotmat';
% rotmat * principle_axis_1 
% pause 
 fun = @(angle) - rotate_off_axis_tau(angle, swimming_axis, D_aligned, Results_body(body_ind).principle_axis_1_index);
% [angle,max_tau] = fminbnd(fun,0,2*pi);
guesses = linspace(0,2*pi,10);  angle = NaN(size(guesses));  max_tau = angle;
for nn = 1:length(guesses)
[angle(nn),temp] = fminsearch(fun,guesses(nn));
max_tau(nn) = -temp;
end
[max_tau,ind] = max(max_tau);
angle = angle(ind);

%     angles = linspace(0,2*pi,200);
%      [d,tau_a] = diagonalness(angles, swimming_axis, D_aligned, Results_body(body_ind).principle_axis_1_index);
%      figure(932)
%      yyaxis left
%      plot(angles*180/pi,d,'bo-','markerfacecolor','b'); grid on
%      yyaxis right
%           plot(angles*180/pi,tau_a,'ro-','markerfacecolor','r'); 
%           hold on
%           plot(repmat(angle*180/pi,1,2),[min(tau_a) max(tau_a)],'k--','linewidth',2)
%           hold off
%           title(['max at angle = ',num2str(angle*180/pi),'    tau = ',num2str(max_tau)]);
%      drawnow
%      
%      pause
     
 
     
 %     Results(results_ind).D_rot = D_rot;
%      Results(results_ind).D_aligned = D_aligned;
%      Results(results_ind).principle_axis_1_index = Results_body(body_ind).principle_axis_1_index;
% %       Results(results_ind).principle_axis_1 = Results_body(body_ind).principle_axis_1;
%         Results(results_ind).rotational_diffusivity = Results_body(body_ind).rotational_diffusivity;
% 
% inds = setdiff(1:3, Results_body(body_ind).principle_axis_1_index);
   Results(results_ind).tau_a_body = max_tau;
%     Results(results_ind).tau_a_body = 1 / (D_aligned(inds(1),inds(1)) + D_aligned(inds(2),inds(2)) );  %timescale for loss of orientation around swimming axis, ignoring fact that this is not the principle axis and there are cross terms.
    % take average value over all possible tensor rotations around swimming
    % axis, since it's unclear how to better deal with this additional issue
%     if isequal(geom_current,[1 0])
%         stopafro
%     end
end