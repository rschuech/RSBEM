V0 = 1;
V = [0.1 0.5 1 5 10];  % um^3, 1 is original used
P0 = 1E-3;  mu = 1E-3 / 1E6;
P = P0 .* V ./ V0;

rows = 150;  cols = 150;
V_mat = repmat(reshape(V,1,1,[]),rows,cols,1);
P_mat = repmat(reshape(P,1,1,[]),rows,cols,1);

sph_rad = (V_mat * 3/4 / pi).^(1/3);

Pareto3.Speed = sqrt( repmat(Pareto_data.Power_eff.Z,1,1,length(V)) .* P_mat ./ sph_rad ./ mu / 6 / pi );

% tau scales with rotational friction coeff
% rotational friction coeff scales with length^3
% tau scales with length^3 or, with volume^1
Pareto3.tau_a = repmat(Pareto_data.tau_a.Z,1,1,length(V)) .* V_mat / V0;

Pareto3.Dm = 1/3 * Pareto3.Speed.^2 .* Pareto3.tau_a;

Pareto3.temporal_SN = Pareto3.Speed .* Pareto3.tau_a.^(3/2);
% pole separation scales with length or volume^(1/3)
Pareto3.pole_separation = repmat(Pareto_data.pole_separation.Z,1,1,length(V)) .* V_mat.^(1/3) / V0.^(1/3);
Pareto3.fore_aft_SN = Pareto3.pole_separation .* sqrt( Pareto3.tau_a );

% Construction Ease always normalized to a sphere, since it's independent
% of size by definition (unlike other tasks which are not normalized now)
Pareto3.construction = repmat(Pareto_data.Construction_Ease.Z,1,1,length(V));

Pareto3.maintenance = 1./ V_mat;  %assume maintenance ease inversely proportional to volume

% Pareto3.tumbling = 

% diffusion mass flow scales with length for a sphere
Pareto3.uptake = repmat(Pareto_data.uptake.Z,1,1,length(V)) .* V_mat.^(1/3) ./ V0^(1/3);


%%


[XP,YP,ZP] = ndgrid(unique(X(:)),unique(Y(:)),V);



goals = {'Speed','Dm','temporal_SN','fore_aft_SN','construction','maintenance'};

goals = {'Speed','Dm','construction','maintenance'};


    perfs = [];
    for g  = 1:length(goals)
        perfs(:,g) = Pareto3.(goals{g})(:);
    end
    
     XPc = XP(:);  YPc = YP(:);  ZPc = ZP(:);
    NaNs = any(isnan(perfs),2);
    perfs(  NaNs , :) = [];
    XPc(NaNs) = [];  YPc(NaNs) = [];  ZPc(NaNs) = [];
    
        [ p, optimal_inds] = paretoFront( perfs );
        
         Optimal = NaN(size(XP));
    for oi = 1:length(optimal_inds)
        Optimal( XP == XPc(optimal_inds(oi)) & YP == YPc(optimal_inds(oi)) & ZP == ZPc(optimal_inds(oi)) ) = 1;
    end
    
    
    