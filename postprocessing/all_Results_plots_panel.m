
load('E:\Hull\Pareto individual tails.mat','Pareto_data');
% load E:\Hull\'Pareto_straight_tail.mat'


% load_measurements;


try, close(104:110), end

shape = 'capsule';
%shape = 'ellipsoid';
shape = 'curved_rod';


AR1_discard = [3.25 3.75 4.25   7.25   1.125 1.25 1.375 1.625 1.75   ...
    2.1 2.2 2.3 2.4 2.6 2.7 2.8 2.9];  %  1.125 1.25 1.375 1.625 1.75
AR2_discard = [  0.12 0.14 0.16 0.18 0.22 0.24 0.26 0.28 0.32 0.34 0.36 0.38 ...
    0.42 0.44 0.46 0.48 0.52 0.54 0.56 0.58 0.62 0.64 0.66 0.68 ...
    0.72 0.74 0.76 0.78 0.82 0.84 0.88 0.94 ...
    ];
%     0.375  ];


discard_particulars = [NaN NaN];

discard_particulars = [discard_particulars;   [ 1.75 0.075; 1.75 0.125; 1.975 0.625; 2.5 0.425; 2.5 0.475; 2.5 0.525; 2.5 0.575; 2.5 0.625; 2.5 0.675; 2.5 0.725; 2.5 0.775; ...
    3 0.825; 3.5 0.175; 3.5 0.425; 3.5 0.475; 3.5 0.525; 3.5 0.575; 3.5 0.625; 3.5 0.675; 3.5 0.725; 3.5 0.775; 3.5 0.825; 3.5 0.8775; ...
    4 0.575; 4 0.625; 4 0.86; 4 0.96; 4.5 0.275; 4.5 0.425; 4.5 0.475; 4.5 0.525; 4.5 0.725; 4.5 0.775; 4.5 0.825; 4.5 0.875; ...
    4.5 0.925; 5 0.86; 5 0.92; 5.5 0.275; 5.5 0.425; 5.5 0.475; 5.5 0.525; 5.5 0.725; 5.5 0.775; 5.5 0.825; 5.5 0.875; 5.5 0.925; ...
    6 0.575; 6 0.86; 6 0.92; 6.5 0.075; 6.5 0.325; 6.5 0.425; 6.5 0.475; 6.5 0.525; 6.5 0.575; 6.5 0.725; 6.5 0.775; ...
    6.5 0.875; 7 0.575;  7 0.86; 7.475 0.375; 7.5 0.86; 9.5 0.075; 9.5 0.325; 11.5 0.175; ...
    1.025 0.055; 1.05 0.1; 1.025 0.2; 1.05 0.3; 2.3 0.65; 3.5 0.875; 3.575 0.95; 3.5 0.925; 4.5 0.575; 5 0.575; 5.5 0.575; 5 0.725; ...
    4.75 0.1; 5.25 0.1;  4.75 0.2;  5.25 0.2;  4.75 0.3;  5.25 0.3;  ...
    4.75 0.5;  5.25 0.5;  4.75 0.55;  5.25 0.55;  4.75 0.575; 5.25 0.575; ...
    4.75 0.6;  5.25 0.6; 4.5 0.625;  5.5 0.625; 4.5 0.675;  5.5 0.675;  4.75 0.7;  5.25 0.7; ...
    4.75 0.75; 5.25 0.75; 4.75 0.8;  5.25 0.8;  4.75 0.86; 5.25 0.86; 4.75 0.9; 5.25 0.9; ...
    4.75 0.92; 5.25 0.92; 4.75 0.96; 5.25 0.96; 4.75 0.98; 6.5 0.925;  4.75 0.4; 5.25 0.4; ...
    1.025 0.05; 1.025 0.1; 1.025 0.15; 1.025 0.25; 1.025 0.3;   1.275 0.4; ...
    1.425 0.45; 1.575 0.5; 1.75 0.55;  2.05 0.65; 2.2 0.7; 2.375 0.75; 2.525 0.8; 2.7 0.85; 2.975 0.9; 1.01 0.2;];   ...
    4.75 0.625; 5 0.625; 5.25 0.625; 4.75 0.65; 5.25 0.65; 4.75 0.675; 5 0.675; 5.25 0.675; ...
    6.5 0.675; 7 0.675; 7.5 0.675; 8 0.675; 8.5 0.675; 9 0.675; 9.5 0.675; 10 0.675; ...
    6.5 0.625; 7 0.625; 7.5 0.625; 8 0.625; 8.5 0.625; 9 0.625; 9.5 0.625; 10 0.625; ];
%     5 0.3; 6.5 0.4; 9 0.85;  ];
%     1.1 0.05; 1.2 0.05; 1.3 0.05; 1.4 0.05; 1.6 0.05; 1.7 0.05; 1.8 0.05; 1.9 0.05; ...
%     1.1 0.15; 1.2 0.15; 1.3 0.15; 1.4 0.185; 1.6 0.15; 1.7 0.15; 1.8 0.15; 1.9 0.15; ...
%     1.1 0.2; 1.2 0.2; 1.3 0.2; 1.4 0.2; 1.6 0.2; 1.7 0.2; 1.8 0.2; 1.9 0.2; ...
%     1.1 0.25; 1.2 0.25; 1.3 0.25; 1.4 0.25; 1.6 0.25; 1.7 0.25; 1.8 0.25; 1.9 0.25; ...
%     1.1 0.3; 1.2 0.3; 1.3 0.3; 1.4 0.3; 1.6 0.3; 1.7 0.3; 1.8 0.3; 1.9 0.3; ...]% ...
%     1.1 0.35; 1.2 0.35; 1.3 0.35; 1.4 0.35; 1.6 0.35; 1.7 0.35; 1.8 0.35; 1.9 0.35; ...
%                         1.3 0.4; 1.4 0.4; 1.6 0.4; 1.7 0.4; 1.8 0.4; 1.9 0.4; ...
%                                           1.6 0.5; 1.7 0.5; 1.8 0.5; 1.9 0.5; ...
%                                                             1.8 0.55; 1.9 0.55; ...
%                                                                       1.9 0.6];

[AR1, AR2] = ndgrid( 1:0.5:12, 0:0.05:1);
AR1_AR2 = roundn([AR1(:) AR2(:)],-10);
skips = ismember(AR1_AR2(:,1), roundn(0.5:1:12,-10)) | ismember(AR1_AR2(:,2), roundn(0.05:0.1:1,-10));
AR1_AR2( ~skips,:) = [];
AR1_AR2 = setdiff(roundn(AR1_AR2,-10), roundn([5 0.65],-10),'rows');  %add global best back in

% 1.125 1.25 1.375 1.625 1.75
straight_discards = [ [1.1:0.1:1.9 2.1:0.1:2.9]' ];  straight_discards = [straight_discards zeros(length(straight_discards),1)];
% AR1_AR2 = [];
more_discards = [3 0.9; 2 0.6;  4.5 0.675;  4.5 0.625;  5.5 0.675;  5.5 0.625;  5 0.625;  5 0.675;  ];



%discard_particulars = [discard_particulars; AR1_AR2; straight_discards;  more_discards];


%%
folder = 'E:\Hull\Results\best\';
files = dir([folder,'Results_all.mat']);

folder = 'E:\Hull\Results\best\';
files = dir([folder,'Results_all_fixed.mat']);  % orig I've been using

% folder = 'E:\Hull\Results\best\';
% files = dir([folder,'Results_all - cleaner.mat']);


% folder = 'E:\Hull\Results\current\';
% files = dir([folder,'Results_straight_tail_all.mat']);

% files = dir([folder,'Results_current.mat']);
% files = dir([folder,'Results_everything_sph_update.mat']);

%  files = dir([folder,'Results_curved_tail_all.mat']);


files = {files.name};
Results = [];
for f = 1:length(files)
    f/length(files)
    file = files{f};
    temp = load([folder,file]);
    try
        results = temp.Results;
    catch
        results = temp.Results0;
    end
    
    
    
    for i = 1:length(results)
        if isempty(results(i).Power_eff)
            results(i).Power_eff = NaN;
        end
    end
    
    
    Results = [Results results];
    
end







% make sure temporal S/N error is fixed
for results_ind = 1:length(Results)
    Results(results_ind).taxis.temporal.SN = Results(results_ind).Adj_Speed * (Results(results_ind).tau_a)^(3/2);
end


Results(ismember(roundn([Results.AR1],-10) , AR1_discard) | ismember(roundn([Results.AR2],-10), AR2_discard)) = [];

Results(ismember( roundn([  [Results.AR1]' [Results.AR2]'  ],-10)  , roundn(discard_particulars, -10)   , 'rows' )) = [];

for i = 1:length(Results)
    Results(i).AR1 = roundn(Results(i).AR1,-10);
    Results(i).AR2 = roundn(Results(i).AR2,-10);
end



%%
[~,inds] = sort([Results.Power_eff],'descend');
Results = Results(inds);
AR1_AR2 = [ [Results.AR1]' [Results.AR2]' ];

AR1s = unique([Results.AR1]);
AR2s = unique([Results.AR2]);
amps = unique([Results.amp]);
lambdas = unique([Results.lambda]);
nlambdas = unique([Results.nlambda]);

%%

taxis = [Results.taxis];
temporal = [taxis.temporal];
temporal_SN = [temporal.SN];
fore_aft = [taxis.fore_aft];
fore_aft_SN = [fore_aft.SN];
buckling = [Results.buckling];
force = [buckling.force];  torque = [buckling.torque];
buckling_force_mean = [force.mean];
buckling_force_max = [force.max];
buckling_torque_mean = [torque.mean];
buckling_torque_max = [torque.max];
% temp = vertcat(Results.fcoeffs);
% temp = vertcat(temp.translation);
% translation = temp(:,1);
% temp = vertcat(Results.fcoeffs);
% temp = vertcat(temp.rotation);
% rotation = temp(:,1);

sph_ind = find( [Results.AR1] == 1 & [Results.AR2] == 0);

for i = 1:length(Results)
    Results(i).temporal_SN = temporal_SN(i);
    
    Results(i).fore_aft_SN = fore_aft_SN(i);
    
    Results(i).lambda_a = Results(i).tau_a * Results(i).Adj_Speed;
    Results(i).buckling_force_mean = buckling_force_mean(i);
    Results(i).buckling_force_max = buckling_force_max(i);
    Results(i).buckling_torque_mean = buckling_torque_mean(i);
    Results(i).buckling_torque_max = buckling_torque_max(i);
    Results(i).tumbling = 1/Results(i).tau_a_body;
    
    
    % recompute tau and everything that depends on it due to fix of method
    swimming_axis = Results(i).path_slope ./ sqrt(sum(Results(i).path_slope.^2));
    
    principle_axis_1 = Results(i).principle_axis_1;
    rotvec = crossprod( principle_axis_1 , swimming_axis);  rotvec = rotvec / sqrt(sum(rotvec.^2));
    angle = acos( dot(swimming_axis, principle_axis_1) );
    rotmat = rotate_arbitrary_vector( rotvec, angle);
    D_aligned = rotmat * Results(i).rotational_diffusivity * rotmat';  % that's how you rotate a tensor
    
    
    
    fun = @(angle) - rotate_off_axis_tau(angle, swimming_axis, D_aligned);
    
    guesses = linspace(0,2*pi,10);  angle = NaN(size(guesses));  max_tau = angle;  % 10 evenly spaced initial guesses seems quite enough (3 would probably do)
    for nn = 1:length(guesses)
        [angle(nn),temp] = fminsearch(fun,guesses(nn));
        max_tau(nn) = -temp;
    end
    [max_tau,ind] = max(max_tau);
    angle = angle(ind);
    
    Results(i).tau_a = max_tau;
    
    
    
    Results(i).Dm = 1/3 * Results(i).Adj_Speed.^2 * Results(i).tau_a;
    
    Results(i).temporal_SN = Results(i).Adj_Speed * (Results(i).tau_a)^(3/2);
    
    Results(i).fore_aft_SN  = Results(i).pole_separation *  sqrt(Results(i).tau_a);
    
    
    %     Results(i).f_ratio = Results(i).fcoeffs.rotation(1) ^ 0.41 / (Results(i).fcoeffs.translation(1) ^2) ;
    % %     Results(i).f_ratio = Results(i).fcoeffs.rotation(1)  - (Results(i).fcoeffs.translation(1) ) ;
    %       tran = (1/Results(i).fcoeffs.translation(1) - min(1./translation) ) / (max(1./translation)-min(1./translation)) ;
    %       rot = (Results(i).fcoeffs.rotation(1) - min(rotation) ) / (max(rotation)-min(rotation)) ;
    %       Results(i).translation_ease = 1 / Results(i).fcoeffs.translation(1) ;
    %       Results(i).rotation_resistance = Results(i).fcoeffs.rotation(1) ;
    
    %       Results(i).f_ratio =  (rot^0 / tran^1) ;
    %       Results(i).f_ratio = Results(i).fcoeffs.rotation(1);
    
end
%%



% titles = {'Swimming Efficiency','Dispersal','Temporal SN','Temporal Chemotaxis','Fore-Aft SN','Fore-Aft Chemotaxis','Construction Ease','Nutrient Uptake'};
% Depvars = {'Power_eff','Dm', 'temporal_SN','temporal_ability','fore_aft_SN','fore_aft_ability','construction'}; %,'tau_a','lambda_a','error_angle'};
Depvars = {'Power_eff','Dm', 'temporal_SN','fore_aft_SN','construction','tumbling','uptake'}; %,'tau_a','lambda_a','error_angle'};
titles = {'Swimming Efficiency','Dispersal','Temporal S/N','Fore-Aft S/N','Construction Ease','Tumbling Ease','Nutrient Uptake'};
names = {'Swimming Efficiency','Dispersal','Temporal SN','Fore-Aft SN','Construction Ease','Tumbling Ease','Nutrient Uptake'};

% Depvars = {'buckling_force_mean','buckling_force_max','buckling_torque_mean','buckling_torque_max'};
% titles = {'Flagellar Buckling Force (mean)','Flagellar Buckling Force (max)','Flagellar Buckling Torque (mean)','Flagellar Buckling Torque (max)'};
% names = {'Flagellar Buckling Force mean','Flagellar Buckling Force max','Flagellar Buckling Torque mean','Flagellar Buckling Torque max'};


% Depvars = {'translation_ease','rotation_resistance'};
% titles = {'translation_ease','rotation_resistance'};
% names = {'translation_ease','rotation_resistance'};
%
% Depvars = {'rotation_resistance'};
% titles = {'rotation_resistance'};
% names = {'rotation_resistance'};

rows = [1 1 2 2 3 3 4 4];
cols = [1 2 1 2 1 2 1 2];

fontsize = 22;
fontsize = 12;

figure(666)
clf
p = panel();
p.marginright = 10;
p.pack(4, 2);
%     p.marginright = 10;
p.de.marginbottom = 8;

dv_end = (length(Depvars) + 1);
for dv = 1:dv_end
    
    
    
    
    
    %     p(rows(dv), cols(dv)).select();
    %  p(rows(dv), cols(dv)).marginright = 20;
    %  p(rows(dv), cols(dv)).marginleft = 20;
    p(rows(dv), cols(dv)).pack('h',{96 []});
    axish = p(rows(dv), cols(dv),1).select();
    p(rows(dv), cols(dv),1).marginright = 1;
    p(rows(dv), cols(dv),2).marginleft = 1;
    
    %  p.de.margin = 1;
    %     figure(dv+41)
    %     clf
    
    npts = 1400;
    npts = 150;
    
    %     npts = 400;
    
    plot_datapts = true;
    %  depvar = 'Power_eff';
    %    depvar = 'Dm';
    %    depvar = 'temporal_SN';
    %    depvar = 'temporal_ability';
    %     depvar = 'fore_aft_SN';
    %    depvar = 'fore_aft_ability';
    %    depvar = 'error_angle';
    %   depvar = 'tau_a';
    %    depvar = 'lambda_a';
    if dv <= length(Depvars)
        depvar = Depvars{dv};
        if ~strcmp(depvar,'construction') && ~strcmp(depvar,'uptake')
            
            var1 = [Results.AR1];   var2 = [Results.AR2];
            X_Y = [[var1]; [var2] ]';  %two vars to use for X and Y axes
            
            [X_Y_unq,inds,inds2] = unique(X_Y,'rows');
            best_depvars = NaN(size(X_Y_unq,1),1);
            best_inds = best_depvars; n_pts = best_depvars;  best_speeds = best_depvars;  best_freqs = best_depvars;
            best_amps = best_depvars;  best_lambdas = best_depvars;  best_nlambdas = best_depvars;
            n_amps = best_depvars;  n_lambdas = best_depvars;  n_nlambdas = best_depvars;
            
            for i = 1:size(X_Y_unq,1) %each unique body shape
                x_y = X_Y_unq(i,:);
                inds_temp = find(X_Y(:,1) == x_y(1) & X_Y(:,2) == x_y(2)); %all the runs with this body shape
                depvars = [Results(inds_temp).(depvar)];
                best_ind = inds_temp( depvars == max(depvars));  %index of best run within this subset
                if ~isempty(best_ind)
                    best_depvars(i) = Results(best_ind).(depvar);
                    
                    best_inds(i) = best_ind(1);
                    n_pts(i) = length(inds_temp);
                    best_speeds(i) = Results(best_ind).Avg_Speed;
                    best_freqs(i) = Results(best_ind).Avg_Omega;
                    best_amps(i) = Results(best_ind).amp;
                    best_lambdas(i) = Results(best_ind).lambda;
                    best_nlambdas(i) = Results(best_ind).nlambda;
                    n_amps(i) = length(unique([Results(inds_temp).amp]));
                    n_lambdas(i) = length(unique([Results(inds_temp).lambda]));
                    n_nlambdas(i) = length(unique([Results(inds_temp).nlambda]));
                else
                    best_depvars(i) = NaN;
                    best_inds(i) = NaN;
                    n_pts(i) = length(inds_temp);
                    best_speeds(i) = NaN;
                    best_freqs(i) = NaN;
                    best_amps(i) = NaN;
                    best_lambdas(i) = NaN;
                    best_nlambdas(i) = NaN;
                    n_amps(i) = NaN;
                    n_lambdas(i)  = NaN;
                    n_nlambdas(i) = NaN;
                end
                
            end
            
        elseif strcmp(depvar,'construction')
            best_depvars =  Pareto_data.Construction_Ease.F(standardize(X_Y_unq, limits));
            % sphere construction ease seems to be NaN so fix that
            %         [~,temp] = ismember(X_Y_unq,[1 0],'rows');  ind = find(temp);
            %         best_depvars(ind) = 1;  %shouldn't need anymore, hopefully fixed
            %         the interpolation (problem was due to including NaNs)
        elseif strcmp(depvar,'uptake')
            best_depvars =  Pareto_data.uptake.F(standardize(X_Y_unq, limits));
            %         % sphere construction ease seems to be NaN so fix that
            %         [~,temp] = ismember(X_Y_unq,[1 0],'rows');  ind = find(temp);
            %         best_depvars(ind) = 1;
        end
        
        
        
        [X,Y] = meshgrid(linspace(min(var1),max(var1),npts),linspace(min(var2),max(var2),npts));
        
        
        %         dep_var = best_depvars / max([Results.(depvar)]);  %normalize to global best performance
        
        [~,temp] = ismember(X_Y_unq,[1 0],'rows');  ind = find(temp);  % ind for sphere, probably first one...
        
        sphere_depvar = best_depvars(ind);
        
        dep_var = best_depvars / sphere_depvar;  %normalize to sphere performance
        
        if strcmp(depvar,'Power_eff')
            power_eff_depvar = dep_var;
        end
        
        % dep_var = log10(max_der');
        %           dep_var =  best_speeds;
        %           dep_var = best_freqs;
        % dep_var = best_depvars;
        %                         dep_var = best_amps;
        %                     dep_var = best_lambdas;
        %          dep_var = best_nlambdas;
        %     dep_var = best_amps ./ best_lambdas;
        %   dep_var = n_pts;
        % dep_var = best_speeds ./ best_freqs;
        %
        
        subset = n_pts >= 1;
        %
        %         method = 'nearest';
        method = 'linear';
        method = 'natural';
        
        [AR1_AR2,inds] = unique([ [Results.AR1]' [Results.AR2]' ],'rows');
        
        AR1_range = max(AR1_AR2(:,1)) - min(AR1_AR2(:,1));
        AR2_range = max(AR1_AR2(:,2)) - min(AR1_AR2(:,2));
        %     limits = [min(AR1_AR2(:,1)) max(AR1_AR2(:,1)); min(AR1_AR2(:,2)) max(AR1_AR2(:,2))];
        limits = [1 10; 0 1];
        %         standardized_X_subset = (  X_Y_unq(subset,1) - min(AR1_AR2(:,1))  ) / AR1_range;
        %         standardized_Y_subset = (  X_Y_unq(subset,2) - min(AR1_AR2(:,2))  ) / AR2_range;
        [standardized_subset] = standardize(X_Y_unq(subset,:), limits);
        
        F = scatteredInterpolant(standardized_subset(:,1), standardized_subset(:,2), dep_var(subset),method,'none');
        %         F = scatteredInterpolant(standardized_X_subset, standardized_Y_subset, dep_var(subset),method,'linear');
        
        
        %         standardized_X = (  X - min(AR1_AR2(:,1))  ) / AR1_range;
        %         standardized_Y = (  Y - min(AR1_AR2(:,2))  ) / AR2_range;
        [standardized_XY] = standardize([X(:) Y(:)], limits);
        standardized_X = NaN(size(X));  standardized_X(:) = standardized_XY(:,1);
        standardized_Y = NaN(size(Y));  standardized_Y(:) = standardized_XY(:,2);
        
        Z = F(standardized_X,standardized_Y);
        
        if strcmp(depvar, 'construction')
            Z = Pareto_data.Construction_Ease.Z;
            Z(isnan(Pareto_data.Power_eff.Z)) = NaN;
        end
        
        %        if strcmp(depvar, 'uptake')
        %         Z = Pareto_data.uptake.Z;
        %         Z(isnan(Pareto_data.Power_eff.Z)) = NaN;
        %     end
        
        if ~strcmp(depvar,'construction') && ~strcmp(depvar,'uptake')
            Pareto_data.(Depvars{dv}).Z = Z;
            Pareto_data.(Depvars{dv}).F = F;
            Pareto_data.(Depvars{dv}).depvar = best_depvars;
        end
        %
        
        fontsizes.labels = 30;
        fontsizes.axes = 22;
        fontsizes.title = 31;
        fontsizes.caxis = 26;
        
        fontsizes.labels = 12;
        fontsizes.axes = 12;
        fontsizes.title = 12;
        fontsizes.caxis = 12;
        
        %the switcheroo is used here so that standardized locations were used to do
        %interpolation, but orig locations are used for plot
        pch = pcolor(X,Y,Z);
        shading interp
        
        set(gca,'fontsize',fontsizes.axes);
        
        if dv == 7
            xlabel('AR_1 (elongation)','fontsize',fontsizes.labels)
            ylabel('AR_2 (curvature)','fontsize',fontsizes.labels)
        elseif ismember(dv,[1 3 5 ])
            set(gca,'XTickLabel',{});
        else
            set(gca,'XTickLabel',{});
            set(gca,'YTickLabel',{});
        end
        
        %
        
        %             xlabel('tail amplitude (\mum)','fontsize',fontsize)
        %         ylabel('tail wavelength (\mum)','fontsize',fontsize)
        %           ylabel('# tail wavelengths','fontsize',fontsize)
        %zlabel('power eff')
        hold on
        
        % plot single best pt
        xtemp = X_Y_unq(subset,1); ytemp = X_Y_unq(subset,2); ztemp =  dep_var(subset);
        %         best_pt =   plot(xtemp(ztemp==max(ztemp(:))),ytemp(ztemp==max(ztemp(:))),'kh','markerfacecolor','k','markersize',12);
        best_pt =   [xtemp(ztemp==max(ztemp(:)))  ytemp(ztemp==max(ztemp(:)))];
        
        
        
        if plot_datapts
            
            % d = plot(X_Y_unq(subset,1),X_Y_unq(subset,2), '.','markersize',0.5,'markerfacecolor','none','markeredgecolor','k');
            switch Depvars{dv}
                %             case 'Power_eff'
                %                 starh = plot(6,0.65,'kp','markersize',20,'markerfacecolor','k');
                %             case 'construction'
                %                 starh = plot(1,0,'kp','markersize',20,'markerfacecolor','k');
                case {'Power_eff','construction','tumbling'}
                    if ~isempty(best_pt)
                        starh = plot(best_pt(1),best_pt(2),'kp','markersize',8,'markerfacecolor','k');
                    end
                otherwise
                    if ~isempty(best_pt)
                        temp = 0.4;
                        starh = plot(best_pt(1),best_pt(2),'p','markersize',8,'markeredgecolor',[temp temp temp],'markerfacecolor',[temp temp temp]);
                    end
                    
            end
        end
        
        hold off
        
        cbar = colorbar(axish);
        
        
        % cbl = cblabel('swimming efficiency / best efficiency','fontsize',fontsize);
        
        % l = legend([d1 d2],'tail crudely optimized','tail not optimized');
        % set(l,'position',[0.58521      0.90343      0.23778     0.067114]);
        %         title(depvar,'interpreter','none')
        title(titles{dv},'fontsize',fontsizes.title,'interpreter','none');
        % cbl = cblabel('swimming speed (\mum/s)','fontsize',fontsize);
        % cbl = cblabel('motor frequency (rev/s)','fontsize',fontsize);
        %  cbl = cblabel('# tail wavelengths','fontsize',fontsize);
        %cbl = cblabel('tail amplitude (\mum)','fontsize',fontsize);
        %    cbl = cblabel('tail wavelength (\mum)','fontsize',fontsize);
        
        
        
        
        
        hold on
        %         pl = plot(best_line.AR1,best_line.AR2,'k-','linewidth',3);
        
        %individually opt tails
        clims = [   0.7  1.042;          1  4.2894;         1  9;            0.0  6.4;          0.6  1;             0.05  1;    0.9 1.4;                     0.9 1.4;                                                       0.95 12;             0.95 12;   1 1.6];
        cvals = {[0.75 0.9 1  1.03],[1.002 1.025 1.05 1.07 1.1 1.5 1.75 2 2.25 2.5 2.7 2.9 3.1 3.3 3.5 3.7 3.9],[1.05 1.5 2.5 4 6 8],[0.1 0.4 0.8 1 1.5 2 3 4 5 6],[0.65 0.7 0.75 0.8 0.9 0.95],[0.1:0.1:1 0.45],[ 0.95 0.975 1 1.05 1.125 1.2 1.25],     [0.95 1 1.1 1.2 1.3 1.4 1.5],    [1.1 1.5 3 2:2:10],    [1.1 1.5 3 2:2:10]  ,[1.1 1.2 1.3 1.4 1.5 1.6] };
        
        %     clims = [ 0.65 1.05];
        %     clims = [ 1 8];
        %     cvals = {linspace(clims(1),clims(2),70)};
        
        %       clims = [ 0.1375  .1525;       1  4.2894;         1  2.0449;          0.0  6.4;          0.6  1;             0.05  1;    0.9 1.4;                     0.9 1.4;                                                       0.95 12;             0.95 12;];
        %     cvals = {[0.75 0.9 1  1.03],[1.1 1.5 2 2.5 3.5],[1.05 1.2 1.5 1.8 2],[0.1 0.4 0.8 1 1.5 2 3 4 5 6],[0.65 0.7 0.75 0.8 0.9 0.95],[0.1:0.1:1],[ 0.95 0.975 1 1.05 1.125 1.2 1.25],     [0.95 1 1.1 1.2 1.3 1.4 1.5],    [1.1 1.5 3 2:2:10],    [1.1 1.5 3 2:2:10] };
        
        
        %straight tail
        %     clims = [ 0.7  1.0423;       1  4.2894;         1  5;          0.0  6.4;          0.6  1;             0.05  1;    0.9 1.4;                     0.9 1.4;                                                       0.95 12;             0.95 12;            1 1.6];
        %     cvals = {[0.75 0.85 0.9 0.95 1  1.015],[1.1 1.2 1.5 2 2.5 3.5],[1.05 1.5 3 4 5],[0.1 0.4 0.8 1 1.5 2 3 4 5 6],[0.65 0.7 0.75 0.8 0.9 0.95],[0.1:0.1:1],[ 0.95 0.975 1 1.05 1.125 1.16 1.2],     [0.95 1 1.1 1.2 1.3 1.4 1.5],    [1.1 1.5 3 5 2:2:10],    [1.1 1.5 3 5 2:2:10]  ,   [1.1 1.2 1.3 1.4 1.5] };
        
        % curved tail
        %     clims = [ 1  1.9;       1  4.2894;         1  8;          0.0  6.4;          0.6  1;             0.05  1;    0.9 1.4;                     0.9 1.4;                                                       0.95 12;             0.95 12;];
        %     cvals = {[1 1.3 1.7 1.8 1.85],[1.1 1.5 2 2.5 3.5],[1.05 1.2 3 4 5 6 7 8 ],[0.1 0.4 0.8 1 1.5 2 3 4 5 6],[0.65 0.7 0.75 0.8 0.9 0.95],[0.1:0.1:1],[ 0.95 0.975 1 1.05 1.125 1.2 1.25],     [0.95 1 1.1 1.2 1.3 1.4 1.5],    [1.1 1.5 3 2:2:10],    [1.1 1.5 3 2:2:10] };
        
        
        % hardcode clims here
        %caxis(clims(dv,:));
        
        %   set(cbar,'fontsize',fontsizes.caxis);  % gets overidden by
        %   cblabel, apparently no way to fix
        %cbl = cblabel('performance relative to spherical body','fontsize',fontsizes.labels);
        
        hold on
        [C,ch] = contour(unique(X),unique(Y),Z,cvals{dv},'k--','linewidth',1);
        
        % clabel(C,ch,'fontsize',fontsizes.axes,'LabelSpacing',1000);
        hold off
        % set(ch,'linewidth',2);
        set(gcf,'Position',[ 680    33   640   976]);
        ylim([0 1]);
        set(gca,'XTick',[1:10])
        p.fontsize = fontsize;
        hold on
        % plot measurements
        %     try, delete(hmeas); end;
        if dv < dv_end
            temp = vertcat(plot_data.mean_unweighted);
            %      hmeas(i) = plot((plot_data(i).AR1),(plot_data(i).AR2),'bo','markerfacecolor','b','markersize',1);
          %  color = [0.8 0.8 0.8];
            hmeas = plot(temp(:,1),temp(:,2),'o','markerfacecolor','k','markeredgecolor','k','markersize',3);
            %         pause
        end
        
        
        % meas = [ [data.AR1]' [data.AR2]'  ];  % all inidividual data points
        % hmeas = plot(meas(:,1),meas(:,2),'ko','markerfacecolor','k','markersize',5);
        % hull = convhull(meas(:,1),meas(:,2));
        % hh = plot(meas(hull,1),meas(hull,2),'k:');
        %
        % edges = {  1:1:26  ,  0:0.05:1   };
        % [N,C] = hist3(meas,edges);
        
        
        hold off
        p(rows(dv), cols(dv),2).select(cbar);
        drawnow
        
        
        
        %               print('-dpng','-r300',['E:\Hull\graphs2\temp2\',names{dv},'  straight tail','.png']);
        % print('-dpng','-r300',['E:\Hull\graphs2\',names{dv},'  curved tail','.png']);
        %       print('-dpng','-r300',['E:\Hull\graphs2\',names{dv},'.png']);
        %           print('-dpng','-r600',['E:\Hull\select_dumps\',titles{dv},'.png']);
        %         pause
        %         saveas(gcf,[outfolder,prefixes{fi},'_',depvar,'_new.png'])
        
        %         saveas(gcf,[outfolder,depvar,'_new.png'])
        
        
        
        
        max_AR2 = 0.9;  %best AR2 not allowed to be above this (to get rid of donuts with long tails being best)
        best_line.AR1 = unique(X_Y_unq(:,1));  %all orig data AR1s
        best_line.AR2 = NaN(size(best_line.AR1));  %best corresponding AR2s to be filled in below
        for i = 1:length(best_line.AR1)
            AR1 = best_line.AR1(i);
            depvar_AR1 = best_depvars(X_Y_unq(:,1) == AR1 & X_Y_unq(:,2) <= max_AR2);
            AR2_AR1 = X_Y_unq(X_Y_unq(:,1) == AR1 & X_Y_unq(:,2) <= max_AR2,2);
            depvars_temp = best_depvars(X_Y_unq(:,1) == AR1 & X_Y_unq(:,2) <= max_AR2);
            %                     [~, max_eff_ind] = max(Z(:,i));
            %                     best_line.AR2(i) = Y(max_eff_ind,1);
            depvar_straight = best_depvars(X_Y_unq(:,1) == AR1 & X_Y_unq(:,2) == 0);
            if numel(AR2_AR1) < 5 || isempty(depvar_straight)
                best_line.AR2(i) = NaN;
                best_line.improvement(i) = NaN;
                if ~isempty(depvar_straight)
                    best_line.straight(i) = depvar_straight;
                    best_line.straight_AR1(i) = AR1;
                else
                    best_line.straight(i) = NaN;
                    best_line.straight_AR1(i) = NaN;
                end
            else
                best_line.AR2(i) = AR2_AR1(depvar_AR1 == max(depvar_AR1));
                best_line.improvement(i) = depvars_temp(depvar_AR1 == max(depvar_AR1))  /  depvar_straight;
                best_line.straight(i) = depvar_straight;
                best_line.straight_AR1(i) = AR1;
            end
        end
        
        best_line.AR1 = best_line.AR1(~isnan(best_line.AR2));
        best_line.improvement = best_line.improvement(~isnan(best_line.AR2));
        
        best_line.AR2 = best_line.AR2(~isnan(best_line.AR2));
        
        
        best_curvatures(dv).AR1 = best_line.AR1;  best_curvatures(dv).AR2 = best_line.AR2;  best_curvatures(dv).improvement = best_line.improvement;  best_curvatures(dv).straight = best_line.straight;   best_curvatures(dv).straight_AR1 = best_line.straight_AR1;
        best_curvatures(dv).AR1 = [1; best_curvatures(dv).AR1];
        best_curvatures(dv).AR2 = [0; best_curvatures(dv).AR2];
        best_curvatures(dv).improvement = [1 best_curvatures(dv).improvement];
        
        %     pause
        
    else
        
        % Last panel is the parameter space plot of shapes, including the optimal shapes
        bodies = [1 0; 10 0; 10 0.95; 1.125 0.35; 2 0.635; 3.5 0.945; 5.5 0.985; 5.5 0; 10 0.65; 6 0.65; 3.5 0.25; 5.5 0.15;   10 0.4;  6 0.35;  3 0.5; 10 0.1; ];
        optimals1 = [1 0; 6 0.65;];
        optimals2 = [10 0; 10 0.1; 10 0.4;];
        for b = 1:size(bodies,1)
            if ismember(bodies(b,:),optimals1,'rows')
                plot(bodies(b,1),bodies(b,2),'kp','markerfacecolor','k','markersize',8);
            elseif ismember(bodies(b,:),optimals2,'rows')
                temp = 0.4;
                plot(bodies(b,1),bodies(b,2),'kp','markersize',8,'markeredgecolor',[temp temp temp],'markerfacecolor',[temp temp temp]);
            else
                plot(bodies(b,1),bodies(b,2),'k*','markerfacecolor','k','markersize',8);
            end
            hold on
        end
        bound_species = plot(observed.species(:,1),observed.species(:,2),'k:','linewidth',1); %species hull
        %         pause
        color  = [0.75 0.75 0.75];
        d = plot(X_Y_unq(subset,1),X_Y_unq(subset,2), '.','markersize',0.5,'markerfacecolor',color,'markeredgecolor',color);
            
        
        
        grid on
        set(gca,'fontsize',fontsizes.axes)
        
        %         xlabel('AR_1 (elongation)','fontsize',fontsizes.labels);
        %         ylabel('AR_2 (curvature)','fontsize',fontsizes.labels);
        title('Parameter Space','fontsize',fontsizes.title);
        
        xlim([1 10]);  ylim([0 1]);
        set(gca,'XTick',[1:10]);
        p.fontsize = fontsize;
        %   set(gca,'XTickLabel',{});
        set(gca,'YTickLabel',{});
        set(gca,'XTick',[1:10]);
        
        
        ax = gca;
        
        dumpfolder = 'E:\Hull\select_dumps\';
        clear axs
        for b = 1:size(bodies,1)
            dumps = dir([dumpfolder,'*AR1_',num2str(bodies(b,1)),'_AR2_',num2str(bodies(b,2)),'_*.mat']);
            dump = dumps.name;
            load([dumpfolder,dump],'Mesh');
            
            [Xf, Yf] = ds2nfu(ax,bodies(b,1) + 0.25, bodies(b,2) - 0.1) ;
            ax_shape = axes('color','none','visible','off','clipping','off','position',[Xf Yf 1/22 1/22]);
            %       set(ax_shape,'visible','on','box','on','xtick','','ytick','');
            axs(b) = ax_shape;
            
            
            orig_limits = [1 10; 0 1];
            shape_limits = [-20 20; -20 20];
            orig2shape.x = polyfit(orig_limits(1,:),shape_limits(1,:),1);
            orig2shape.y = polyfit(orig_limits(2,:),shape_limits(2,:),1);
            shape_coord = [polyval(orig2shape.x,bodies(b,1)) polyval(orig2shape.y,bodies(b,2)) 0]';
            %      Mesh = shiftMesh(Mesh,[shape_coord - Mesh(2).refpoints(:,1)]);
            
            %        figure(453)
            [s,e] = plot_mesh(Mesh);
            set(e,'edgealpha',0.0);
            xlim([Mesh(2).refpoints(1,1) Mesh(2).refpoints(1,1) + 4])
            ylim([Mesh(2).refpoints(2,1) Mesh(2).refpoints(2,1) + 2])
            zlim([-1 3.5])
            
            % hold on
            light
            drawnow
            % pause
        end
        
        %     axes(ax);
        
        %         for i = 1:length(plot_data)
        %             %      hmeas(i) = plot((plot_data(i).AR1),(plot_data(i).AR2),'bo','markerfacecolor','b','markersize',1);
        %             hmeas(i) = plot((plot_data(i).mean_unweighted(1)),(plot_data(i).mean_unweighted(2)),'ko','markerfacecolor','k','markersize',2);
        %         end
        
    end
    
    
end  %depvars



%%

return






%     best_curvatures(end+1).AR1 = best_curvatures(1).AR1;   best_curvatures(end).AR2 = best_curvatures(1).AR2;

%%
best_curvatures(1).options.AR2 = slmset('knots',12,'decreasing',[3.5 10]);  best_curvatures(1).options.improvement = slmset('knots',3);
best_curvatures(2).options.AR2 = slmset('knots',6,'concaveup',[3 10]);  best_curvatures(2).options.improvement = slmset('knots',3);
best_curvatures(4).options.AR2 = slmset('knots',6,'concaveup',[3 10]);  best_curvatures(4).options.improvement = slmset('knots',3);
best_curvatures(6).options.AR2 = slmset('knots',3,'increasing',[2 10]);  best_curvatures(6).options.improvement = slmset('knots',3);



%%

fontsizes.labels = 26;
fontsizes.axes = 18;
fontsizes.title = 29;
fontsizes.caxis = 26;
fontsizes.legend = 22;



figure(987)
clf
styles = {'--',':','-^','-v','-.'};
ii = 0;
for i = [1 2 3 4 6]  %skip SN
    ii = ii + 1;
    yyaxis left
    
    %         slm = slmengine(best_curvatures(i).AR1,best_curvatures(i).AR2, best_curvatures(i).options.AR2);
    %         AR1_refined = linspace(1,10,200);
    %         splined = slmeval(AR1_refined,slm);
    %         plot(AR1_refined,splined,styles{ii},'linewidth',2);  hold on;
    plot(best_curvatures(i).AR1,best_curvatures(i).AR2,styles{ii},'linewidth',2);
    hold on
    yyaxis right
    
    %         slm = slmengine(best_curvatures(i).AR1,best_curvatures(i).improvement, best_curvatures(i).options.improvement);
    %         AR1_refined = linspace(1,10,200);
    %         splined = slmeval(AR1_refined,slm);
    %         plot(AR1_refined,(splined - 1)*100,styles{ii},'linewidth',2);
    plot(best_curvatures(i).AR1,(best_curvatures(i).improvement - 1)*100,styles{ii},'linewidth',2);
    hold on
    drawnow
    %          pause
end
yyaxis left
ihans = plot([1],[0],'k--',[1],[0],'k:',[1],[0],'k-^',[1],[0],'k-v',[1],[0],'k-.');
set(ihans,'linewidth',2.5);

set(gca,'fontsize' ,fontsizes.axes);

%    set(ihans,'Visible','off');
%  leg =   legend(ihans,{'Swimming Efficiency','Dispersal','Temporal Chemotaxis','Fore-Aft Chemotaxis'},'location','best');
leg = legendflex(ihans,{'Swimming Efficiency','Dispersal','Temporal S/N','Fore-Aft S/N','Tumbling Ease'} , 'nrow',2,'fontsize',fontsizes.legend,'xscale',1.5);
% set(leg,'fontsize',fontsizes.labels);
%  set(leg,'position',[  0.3924        0.747      0.27969      0.17247]);
set(leg,'position',[       144.25        752.8       672.75       82.027]);
hold off

title('Benefits of Curvature','fontsize',fontsizes.title);


yyaxis left
ylim([0 0.9])
ylabel('best AR_2 (curvature)','fontsize',fontsizes.labels);
yyaxis right
ylim([0 Inf])
ylabel('% improvement over straight','fontsize',fontsizes.labels)
xlabel('AR_1 (elongation)','fontsize',fontsizes.labels)
set(gca,'XTick',1:10);
grid on

print('-dpng','-r300',['E:\Hull\graphs2\','curvature_improvement  curved tail','.png']);

%%

return
% end  %files
%% add Constrution Ease to pareto_data struct
[mean_curv, max_curv, mean_abs_curv, mean_max_curv] = curved_rod_curvature(X, Y, 1);
Z = 1./mean_abs_curv;  Z = Z/max(Z(:));
F = scatteredInterpolant(standardized_X(~isnan(Z)), standardized_Y(~isnan(Z)), Z(~isnan(Z)),method,'none');
Pareto_data.Construction_Ease.Z = Z;   Pareto_data.Construction_Ease.F = F;
[mean_curv2, max_curv2, mean_abs_curv2, mean_max_curv2] = curved_rod_curvature(X_Y_unq(:,1), X_Y_unq(:,2), 1);
Z = 1./mean_abs_curv2;  Z = Z/max(Z(:));
Pareto_data.Construction_Ease.depvar = Z;


%%
sweeps_already_done = [eff.curved_rod.AR1s' eff.curved_rod.AR2s' eff.curved_rod.amps' eff.curved_rod.lambdas' eff.curved_rod.nlambdas'];

%% straight rod line graph

straight_inds = find(X_Y_unq(:,2) == 0);
straight_bodies = X_Y_unq(straight_inds,:);

% straight_dep_vars = dep_var(straight_inds);
clear straight_dep_vars
straight_dep_vars.opt = straight_effs.optimized_tail(straight_inds);
straight_dep_vars.straight = straight_effs.straight_tail(straight_inds);
straight_dep_vars.curved = straight_effs.curved_tail(straight_inds);

figure(932)
plot(straight_bodies(:,1),straight_dep_vars.opt,'bo-','markerfacecolor','b'); hold on
plot(straight_bodies(:,1),straight_dep_vars.straight,'rv-','markerfacecolor','r');
plot(straight_bodies(:,1),straight_dep_vars.curved,'r^-','markerfacecolor','r');
grid on
xlabel('AR_1','fontsize',fontsize-3)
ylabel('Normalized Power eff','fontsize',fontsize-3)
title('Straight Rod Swimming Efficiency','fontsize',fontsize)
xlim([1 Inf]);
ylim([0.6 1.35]);
set(gca,'fontsize',fontsize-5)
hold on
plot([1.67 1.67],[0.6 1.35],'--k');
hold off
legend('individually optimized tails','optimal straight rod tail','optimal curved rod tail',...
    'Shum et al ellipsoid optimum','location','best')

% print tail for best straight rod (read AR1 off of straight rod line
% graph)
ind = find(X_Y_unq(:,1)==1.5 & X_Y_unq(:,2)==0)
ind = 15;  [best_amps(ind) best_lambdas(ind)  best_nlambdas(ind)]

%  0.45577       3.2192       1.3104

% print tail for best curved rod (read AR1,AR2 off of main graph)
ind = find(X_Y_unq(:,1)==5 & X_Y_unq(:,2)==0.65)
[best_amps(ind) best_lambdas(ind)  best_nlambdas(ind)]

%   0.61746       4.2255       1.3435

%%
bodies0 = bodies;  guesses0 = guesses;
%%
racks = [1 1 2 2 ];
reps = [1 2 1 2 1 2 1 2];
allbods = [];
for n = 1:length(racks)
    % n = 4;
    bodies = bodies0(n:length(racks):end,:);
    guesses = guesses0(n:length(racks):end,:);
    inds = randperm(size(bodies,1));
    bodies = bodies(inds,:);
    guesses = guesses(inds,:);
    save(['E:\Hull\sweep_CFD0',num2str(racks(n)),'_',num2str(reps(n)),'.mat'],'bodies','guesses');
    allbods = [allbods; bodies];
end
%% new_sweep to run same tail for all bodies
%   0.61746       4.2255       1.3435
clear new_sweep

new_sweep.AR1 = X_Y_unq(:,1);
new_sweep.AR2 = X_Y_unq(:,2);
new_sweep.amp = repmat(0.61746,size(X_Y_unq,1),1);
new_sweep.lambda = repmat( 4.2255 ,size(X_Y_unq,1),1);
new_sweep.nlambda = repmat(1.3435,size(X_Y_unq,1),1);

save('E:\Hull\new_sweep_opt_tails.mat','new_sweep');
