


%%

nSF1 = 300;  nSF2 = 300;

nSF1 = 60;  nSF2 = 40;

nSF1 = 120;  nSF2 = 80;

% nSF1 = 1200;  nSF2 = 800;
nSF1 = 240;  nSF2 = 160;

nSF1 = 125;  nSF2 = 100;  % "final"

nSF1 = 1250;  nSF2 = 1000; % for smoothed headliner


dA = (10 - 1)/(nSF1 - 1) * (1 - 0)/(nSF2 - 1);  %area of a rectangle centered on each interrogation point

xp = linspace(1,10,nSF1);  yp = linspace(0,1,nSF2);
[XP,YP] = ndgrid(linspace(1,10,nSF1),linspace(0,1,nSF2));
[standardized_XPYP] = standardize([XP(:) YP(:)], [1 10; 0 1]);
standardized_XP = NaN(size(XP));  standardized_XP(:) = standardized_XPYP(:,1);
standardized_YP = NaN(size(YP));  standardized_YP(:) = standardized_XPYP(:,2);


% Pareto_data.interped = NaN(nSF1,nSF2,length(Fs));     Pareto_data.rank =  NaN(nSF1,nSF2,length(Fs));
interped_temp = NaN(nSF1,nSF2,length(Fs));  ranked_temp = NaN(nSF1,nSF2,length(Fs));

parfor f = 1:length(Fs)
    interped_temp(:,:,f) = Fs{f}(standardized_XP,standardized_YP);
    temp = interped_temp(:,:,f);
    [~,I] = sort( temp(:) ,'ascend' );
    inds = 1:length(I);
    temp = ranked_temp(:,:,f);
    temp(I) = inds;  %I rearranges inds so that we get the sorted ranks
    
    temp(isnan(interped_temp(:,:,f))) = NaN;
    ranked_temp(:,:,f) = temp;
end

Pareto_data.interped = interped_temp;  Pareto_data.rank = ranked_temp;

% run when you need to calc uber refined interped data for coloring smoothed Pareto
% regions
if nSF1 > 500
    Pareto_data_refined = Pareto_data;
    Pareto_data_refined.nSF1 = nSF1;  Pareto_data_refined.nSF2 = nSF2;
    Pareto_data_refined.xp = xp;  Pareto_data_refined.yp = yp;
    Pareto_data_refined.XP = XP;  Pareto_data_refined.YP = YP;
    return
end
% only need to run till here for smoothed Pareto data - no need to do
% Pareto calc below!
stopa
%%

% goal_inds = setdiff(1:length(Pareto_data.fullnames) , [1 3]);
% goal_inds = setdiff(1:length(Pareto_data.fullnames) , []);

% goal_inds = [1 3 5 6 7:12  18:23 30 ];

goal_inds = [1 3 5 6 9 20 22]; % only three logically reasonable construction eases

% construction_inds = [7:12  18:23 30];  % no longer considering surface max based variants
construction_inds = [9 20 22 ]; 

goal_combos = {};
% loop = length(goal_inds):-1:3;
% loop = 4:-1:3;
% loop = [  6 5 4 3 ];
loop = [5 4 3];

for ngoals = loop
    
    
    temp = combnk(goal_inds,ngoals);
    
    
    for t = 1:size(temp,1)
        
        % don't include more than 1 construction ease variant
        ism = ismember(construction_inds, temp(t,:));
        if sum(ism) > 1
            continue
        end
        
        if ismember( temp(t,:)  , [2 4]) % fore-aft, Dm
%             continue
        end
        
        goal_combos = [goal_combos; {temp(t,:)}];
        
    end
end

goal_combos = flipud(goal_combos);

%  goal_combos = {[5 6 18]};
% goal_combos = {[5 6 19]};
% goal_combos = { [5 6] ,[5 22] ,[6 22]};

% goal_combos = {[5 6 7],[5 6 8],[5 6 9],[5 6 10],[5 6 11],[5 6 12],[5 6 13],[5 6 14],[5 6 15],[5 6 16],[5 6 17],[5 6 18],[5 6 19],[5 6 20]};
% goal_combos = {[5 6 22],... eff, SNR, weighted abs gauss curv
%    ... [5 6 11],... eff, SNR, unweighted abs gauss curv
%     [3 6 22],... tumbling, SNR, weighted abs gauss curv (unweighted is same)
%     [3 5 6 22],... eff, SNR, tumbling, weighted abs gauss curv (unweighted is same)
%     [1 3 22],... uptake, tumbling, weighted abs gauss curv    (adding SNR is same, most/all constructions are same)
%     [1 5 22],... uptake, eff, weighted abs gauss curv (most/all constructions are same)
%     };

goal_combos = {[5 6 22]};
%%
Pareto_precomputed = struct('optimal_inds',cell(length(goal_combos),1),...
    'Optimal',cell(length(goal_combos),1),...
    'A_pareto',cell(length(goal_combos),1),...
    'goals',cell(length(goal_combos),1),...
    'fullnames',cell(length(goal_combos),1));

ppm = ParforProgressStarter2('precomputing Pareto regions', length(goal_combos));


for go = 1:length(goal_combos)
%     go / length(goal_combos)
    
    goals = goal_combos{go};
    %     Pareto_data.fullnames{goals}
    
    
    
    perfs = NaN([size(Pareto_data.interped,1) * size(Pareto_data.interped,2) , length(goals)]);
    for g  = 1:length(goals)
        perfs(:,g) = reshape( Pareto_data.interped(:,:,goals(g)) ,[],1);
    end
    
    XPc = XP(:);  YPc = YP(:);
    NaNs = any(isnan(perfs),2);
    perfs(  NaNs , :) = [];
    XPc(NaNs) = [];  YPc(NaNs) = [];
    
    
    temp = round(perfs,8,'significant');
    [ p, optimal_inds] = paretoFront( temp );
    
    
    
    Optimal = NaN(size(XP));
    for oi = 1:length(optimal_inds)
        Optimal( XP == XPc(optimal_inds(oi)) & YP == YPc(optimal_inds(oi)) ) = 1;
    end
    
    A_pareto = dA * length(optimal_inds);  % approx area of Pareto region
    
    
    Pareto_precomputed(go).optimal_inds = optimal_inds;
    Pareto_precomputed(go).Optimal = Optimal;
    Pareto_precomputed(go).A_pareto = A_pareto;
    Pareto_precomputed(go).goals = goals;
    Pareto_precomputed(go).fullnames = Pareto_data.fullnames(goals);
    
    ppm.increment(go);
    
end

delete(ppm);