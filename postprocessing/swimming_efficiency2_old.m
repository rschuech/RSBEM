

shape = 'capsule';
%shape = 'ellipsoid';
shape = 'curved_rod';


AR1_discard = [3.25 3.75 4.25   7.25       ];  %  1.125 1.25 1.375 1.625 1.75
AR2_discard = [  0.12 0.14 0.16 0.18 0.22 0.24 0.26 0.28 0.32 0.34 0.36 0.38 ...
    0.42 0.44 0.46 0.48 0.52 0.54 0.56 0.58 0.62 0.64 0.66 0.68 ...
    0.72 0.74 0.76 0.78 0.82 0.84 0.88 0.94 ];


discard_particulars = [NaN NaN];

discard_particulars = [discard_particulars;   [1.2 0.375; 1.75 0.075; 1.75 0.125; 1.975 0.625; 2.5 0.425; 2.5 0.475; 2.5 0.525; 2.5 0.575; 2.5 0.625; 2.5 0.675; 2.5 0.725; 2.5 0.775; ...
    3 0.825; 3.5 0.175; 3.5 0.425; 3.5 0.475; 3.5 0.525; 3.5 0.575; 3.5 0.625; 3.5 0.675; 3.5 0.725; 3.5 0.775; 3.5 0.825; 3.5 0.8775; ...
    4 0.575; 4 0.625; 4 0.86; 4 0.96; 4.5 0.275; 4.5 0.425; 4.5 0.475; 4.5 0.525; 4.5 0.725; 4.5 0.775; 4.5 0.825; 4.5 0.875; ...
    4.5 0.925; 5 0.86; 5 0.92; 5.5 0.275; 5.5 0.425; 5.5 0.475; 5.5 0.525; 5.5 0.725; 5.5 0.775; 5.5 0.825; 5.5 0.875; 5.5 0.925; ...
    6 0.575; 6 0.86; 6 0.92; 6.5 0.075; 6.5 0.325; 6.5 0.425; 6.5 0.475; 6.5 0.525; 6.5 0.575; 6.5 0.625; 6.5 0.675; 6.5 0.725; 6.5 0.775; ...
    6.5 0.875; 7 0.575; 7 0.625; 7 0.86; 7.475 0.375; 7.5 0.86; 9.5 0.075; 9.5 0.325; 11.5 0.175; ...
    1.025 0.055; 1.05 0.1; 1.025 0.2; 1.05 0.3; 2.3 0.65; 3.5 0.875; 3.575 0.95; 3.5 0.925; 4.5 0.575; 5 0.575; 5.5 0.575; 5 0.725; ...
      4.75 0.1; 5.25 0.1;  4.75 0.2;  5.25 0.2;  4.75 0.3;  5.25 0.3;  ...
    4.75 0.5;  5.25 0.5;  4.75 0.55;  5.25 0.55;  4.75 0.575; 5.25 0.575; ...
    4.75 0.6;  5.25 0.6; 4.5 0.625;  5.5 0.625; 4.5 0.675;  5.5 0.675;  4.75 0.7;  5.25 0.7; ...
    4.75 0.75; 5.25 0.75; 4.75 0.8;  5.25 0.8;  4.75 0.86; 5.25 0.86; 4.75 0.9; 5.25 0.9; ...
    4.75 0.92; 5.25 0.92; 4.75 0.96; 5.25 0.96; 4.75 0.98; 6.5 0.925;  4.75 0.4; 5.25 0.4; ...
    1.025 0.05; 1.025 0.1; 1.025 0.15; 1.025 0.25; 1.025 0.3;  1.1 0.35; 1.275 0.4; ...
    1.425 0.45; 1.575 0.5; 1.75 0.55; 1.9 0.6; 2.05 0.65; 2.2 0.7; 2.375 0.75; 2.525 0.8; 2.7 0.85; 2.975 0.9; 1.01 0.2;];   ...
    4.75 0.625; 5 0.625; 5.25 0.625; 4.75 0.65; 5.25 0.65; 4.75 0.675; 5 0.675; 5.25 0.675;]% ...
    %[5.5 0.65; 1.5 0.3; 7.5 0.55; 7.5 0.45;  8 0.4; 6.5 0.65; 8.5 0.9; 10 0.5; ];];
    %[7.5 0.55; 8 0.4;  
%     5.75 0.625; 6 0.625; 6.25 0.625; 5.75 0.65; 6.25 0.65; 5.75 0.675; 6 0.675; 6.25 0.675; 12 0.95; ...
%     5 0.35; 3.5 0.4; 11.5 0.65; 12 0.15; 7.5 0.55; 2.5 0.05; 5.5 0.55; 4 0.6; 5.5 0.65; 5.5 0.9; 6 0.9; 7 0.8; 8 0.75; 6.5 0.35; 11 0.6; 11.5 0.2;];

[AR1, AR2] = ndgrid( 1:0.5:12, 0:0.05:1);
AR1_AR2 = roundn([AR1(:) AR2(:)],-10);
skips = ismember(AR1_AR2(:,1), roundn(0.5:1:12,-10)) | ismember(AR1_AR2(:,2), roundn(0.05:0.1:1,-10));
AR1_AR2( ~skips,:) = [];
AR1_AR2 = setdiff(roundn(AR1_AR2,-10), roundn([5 0.65],-10),'rows');  %add global best back in


straight_discards = [ [1.1:0.1:1.9 2.1:0.1:2.9]' ];  straight_discards = [straight_discards zeros(length(straight_discards),1)];
% AR1_AR2 = [];
more_discards = [3 0.9; 2 0.6;  4.5 0.675;  4.5 0.625;  5.5 0.675;  5.5 0.625;  5 0.625;  5 0.675;  ];



%discard_particulars = [discard_particulars; AR1_AR2; straight_discards;  more_discards];



%%
Results = [];

temp = load(['E:\Hull\Results\results0','.mat']);
Results = temp.Results0;


for rack = 1:4
    temp = load(['E:\Hull\Results\results2_CFD0',num2str(rack),'.mat']);
    for i = 1:length(temp.Results)
        temp.Results(i).rack = rack;
    end
    Results = [Results temp.Results];
    
end


Results(ismember(roundn([Results.AR1],-10) , AR1_discard) | ismember(roundn([Results.AR2],-10), AR2_discard)) = [];

Results(ismember( roundn([  [Results.AR1]' [Results.AR2]'  ],-10)  , discard_particulars   , 'rows' )) = [];

for i = 1:length(Results)
    Results(i).AR1 = roundn(Results(i).AR1,-10);
    Results(i).AR2 = roundn(Results(i).AR2,-10);
end


%%
folder = 'E:\Hull\Results\current\';

files = dir([folder,'*.mat']);
files = dir([folder,'Results_all.mat']);

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
    
    if ~isfield(results,'Power_eff')
        for i = 1:length(results)
            results(i).Avg_Power = results(i).Avg_Omega * results(i).Motor_Torque;  %one of these should be constant and the other varying unless I eventually implement a constant power motor condition
            
            mu = 1E-9; %kg / (s micron)
            
            V = 1 ;  %all shapes will have this volume = 1 um^3
            sphererad = (V*3/4/pi)^(1/3); %only used to get tail radius according to Shum et al
            
            results(i).Torque_eff = 8*pi*mu*sphererad^2*results(i).Avg_Speed / results(i).Motor_Torque;
            
            results(i).Power_eff = 6*pi*mu*sphererad*results(i).Avg_Speed^2 / results(i).Avg_Power;
        end
    end
    
    for i = 1:length(results)
        if isempty(results(i).Power_eff)
            results(i).Power_eff = NaN;
        end
    end
    
%     try
%         results = rmfield(results,'fcoeffs');
%         results = rmfield(results,'rack');
%     end
%     try
%         results = rmfield(results,{'Adj_Speed','forced','D'});
%     end
%     try
%         results = rmfield(results,{'tau_a','Dm','error_angle'});
%     end
%     try
%         results = rmfield(results,{'pole_separation','taxis'});
%     end
    
    Results = [Results results];
    
end

Results(ismember(roundn([Results.AR1],-10) , AR1_discard) | ismember(roundn([Results.AR2],-10), AR2_discard)) = [];

Results(ismember( roundn([  [Results.AR1]' [Results.AR2]'  ],-10)  , roundn(discard_particulars, -10)   , 'rows' )) = [];

for i = 1:length(Results)
    Results(i).AR1 = roundn(Results(i).AR1,-10);
    Results(i).AR2 = roundn(Results(i).AR2,-10);
end
%%

for i = 1:length(Results)
    
    %   if ismember(AR1, AR1_discard) || ismember(AR2, AR2_discard)
    %       continue
    %   end
    
    
    Results(i).Avg_Power = Results(i).Avg_Omega * Results(i).Motor_Torque;  %one of these should be constant and the other varying unless I eventually implement a constant power motor condition
    
    mu = 1E-9; %kg / (s micron)
    
    V = 1 ;  %all shapes will have this volume = 1 um^3
    sphererad = (V*3/4/pi)^(1/3); %only used to get tail radius according to Shum et al
    
    Results(i).Torque_eff = 8*pi*mu*sphererad^2*Results(i).Avg_Speed / Results(i).Motor_Torque;
    
    Results(i).Power_eff = 6*pi*mu*sphererad*Results(i).Avg_Speed^2 / Results(i).Avg_Power;
    
    
end

%%
folder = 'E:\Hull\Results\';
outfolder = 'E:\Hull\graphs\';
prefixes = {'straight_opt_tail','curved_opt_tail','best_opt_tail'};
files = {'results_straight_tail_all_variable_interp_pts_CFD04.mat' , 'results_curved_tail_all_4_interp_pts_CFD04.mat' , 'results_opt_tails_all_variable_interp_pts_CFD04.mat'};

for fi = 1:length(files)
    file = files{fi};
    
    load([folder,file]);
    
  
    
    
    %%
    [~,inds] = sort([Results.Power_eff],'descend');
    Results = Results(inds);
    
    
    AR1s = unique([Results.AR1]);
    AR2s = unique([Results.AR2]);
    amps = unique([Results.amp]);
    lambdas = unique([Results.lambda]);
    nlambdas = unique([Results.nlambda]);
    %%
    
    taxis = [Results.taxis];
    temporal = [taxis.temporal];
    temporal_SN = [temporal.SN];
    temporal_ability = [temporal.ability];
    fore_aft = [taxis.fore_aft];
    fore_aft_SN = [fore_aft.SN];
    fore_aft_ability = [fore_aft.ability];
    for i = 1:length(Results)
        Results(i).temporal_SN = temporal_SN(i);
        Results(i).temporal_ability = temporal_ability(i);
        Results(i).fore_aft_SN = fore_aft_SN(i);
        Results(i).fore_aft_ability = fore_aft_ability(i);
        Results(i).lambda_a = Results(i).tau_a * Results(i).Adj_Speed;
    end
   %%
    
    
    
    
    Depvars = {'Power_eff','Dm', 'temporal_SN','temporal_ability','fore_aft_SN','fore_aft_ability','tau_a','lambda_a','error_angle'};
    
    for dv = 1:length(Depvars)
       %%
        depvar = Depvars{dv};
        
        fontsize = 22;
        
        figure(dv+103)
        
        npts = 1400;
        %npts = 20;
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
        
        
        
        var1 = [Results.AR1];   var2 = [Results.AR2];
        % var1 = [Results.amp];  var2 = [Results.nlambda];
        %var1 = eff.curved_rod.nlambdas;  var2 = eff.curved_rod.amps;
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
        
        
        
        
        [X,Y] = meshgrid(linspace(min(var1),max(var1),npts),linspace(min(var2),max(var2),npts));
        dep_var = best_depvars / max([Results.(depvar)]);
        % dep_var = log10(max_der');
        %           dep_var =  best_speeds;
        %           dep_var = best_freqs;
       % dep_var = best_depvars;
%                        dep_var = best_amps;
%                     dep_var = best_lambdas;
%          dep_var = best_nlambdas;
        %   dep_var = n_pts;
        % dep_var = best_speeds ./ best_freqs;
        %

        subset = n_pts >= 1;
        % 
        method = 'nearest';
         method = 'linear';
%               method = 'natural';
                          
        [AR1_AR2,inds] = unique([ [Results.AR1]' [Results.AR2]' ],'rows');
        
        AR1_range = max(AR1_AR2(:,1)) - min(AR1_AR2(:,1));
        AR2_range = max(AR1_AR2(:,2)) - min(AR1_AR2(:,2));
        %  standardized_dists = sqrt(  1/AR1_range^2*( X_Y_unq(subset,1) ).^2  +  1/AR2_range^2*( X_Y_unq(subset,2) ).^2  );
        standardized_X_subset = (  X_Y_unq(subset,1) - min(AR1_AR2(:,1))  ) / AR1_range;
        standardized_Y_subset = (  X_Y_unq(subset,2) - min(AR1_AR2(:,2))  ) / AR2_range;
        
        %         F = scatteredInterpolant(X_Y_unq(subset,1), X_Y_unq(subset,2), dep_var(subset),method,'none');
        
        F = scatteredInterpolant(standardized_X_subset, standardized_Y_subset, dep_var(subset),method,'none');
%         F = scatteredInterpolant(standardized_X_subset, standardized_Y_subset, dep_var(subset),method,'linear');
        %F = scatteredInterpolant(ARs_mirrored',AR2s_mirrored',torque_effs_mirrored'/max(torque_effs_mirrored),'natural','none');
        
        %         Z = F(X,Y);
        standardized_X = (  X - min(AR1_AR2(:,1))  ) / AR1_range;
        standardized_Y = (  Y - min(AR1_AR2(:,2))  ) / AR2_range;
        Z = F(standardized_X,standardized_Y);
        
        
        
        
        %     figure(456)
        %     surf(X,Y,Z)
        %     xlabel('AR')
        %     ylabel('AR2')
        %     zlabel('power eff')
        %     hold on
        %     d = plot3(ARs_mirrored,AR2s_mirrored,power_effs_mirrored / max(power_effs_mirrored),'o','markerfacecolor','k');
        %     %d = plot3(ARs_mirrored,AR2s_mirrored,torque_effs_mirrored / max(torque_effs_mirrored),'o','markerfacecolor','k');
        %
        fontsize = 22;
        %
        %    Loading body mesh     C:\Users\rudi\Desktop\RD\swept_meshes\curved_rod_AR1_11.5_AR2_0.35.dat
        % Loading tail mesh     C:\Users\rudi\Desktop\RD\swept_meshes\tail_radius_0.03101752454497_amp_1.26629518135281_lambda_7.838700515612396_nlambda_2.797322544864104.dat
        
        
        %         pcolor(X,Y,Z);
        %the switcheroo is used here so that standardized locations were used to do
        %interpolation, but orig locations are used for plot
        pcolor(X,Y,Z)
        % surf(X,Y,Z)
        shading interp
        xlabel('AR_1','fontsize',fontsize)
        ylabel('AR_2','fontsize',fontsize)
        %
        %             xlabel('tail amplitude (\mum)','fontsize',fontsize)
        %         ylabel('tail wavelength (\mum)','fontsize',fontsize)
        %           ylabel('# tail wavelengths','fontsize',fontsize)
        %zlabel('power eff')
        hold on
        
        % plot single best pt
        xtemp = X_Y_unq(subset,1); ytemp = X_Y_unq(subset,2); ztemp =  dep_var(subset);
        best_pt =   plot(xtemp(ztemp==max(ztemp(:))),ytemp(ztemp==max(ztemp(:))),'kh','markerfacecolor','k','markersize',12);
        temp = caxis;
        %         c95 = contour(X,Y,Z,[0.95 0.95],'k-');
        caxis(temp);
        
        % caxis([0.6 1]);
%         clear c
%         [~,c(1)] = contour(X,Y,Z,[0.90 0.90],'k:');
%          [~,c(2)] = contour(X,Y,Z,[0.95 0.95],'k--');
%         [~,c(3)] = contour(X,Y,Z,[0.99 0.99],'k-');
        


%  vec = linspace(1,1.5,5);  vals = F((  vec - min(AR1_AR2(:,1))  ) / AR1_range, (  zeros(size(vec)) - min(AR1_AR2(:,2))  ) / AR2_range);
% % 
% for cl = [0.6:0.025:0.95   vals  0.96:0.01:0.99 0.992:0.002:0.999999999999   0.93 0.94 0.972:0.001:0.978 0.982:0.002:0.988  0.999]
%     [temp] = contourcs(unique(X),unique(Y),Z,[cl cl]);
% %     lengths = [temp.Length];
% %     [~,ind] = max(lengths);
% %     if ~isempty(ind)
% for ind = 1:length(temp)
%     plot(temp(ind).X,temp(ind).Y,'k--');
% %    cl
% %    drawnow
% %    pause
%      end
% end
nth = 5;

c_vals = F(standardized_X_subset(1:nth:end), standardized_Y_subset(1:nth:end)); %draw a contour through every data point

for cl = c_vals'
    [temp] = contourcs(unique(X),unique(Y),Z,[cl cl]);
    for ind = 1:length(temp)
        plot(temp(ind).X,temp(ind).Y,'k--');
        %    cl
        %    drawnow
        %    pause
    end
end



        % con = contour(X,Y,Z,[0.98 0.98],'k-');
        % clabel(con);
        if plot_datapts
            try
                fucko
                d1 = plot(X_Y_unq(subset & isnan([next_iter.amp])',1),X_Y_unq(subset & isnan([next_iter.amp])',2), '.','markersize',8,'markerfacecolor','none','markeredgecolor','k');
                d2 = plot(X_Y_unq(subset & ~isnan([next_iter.amp])',1),X_Y_unq(subset & ~isnan([next_iter.amp])',2), 'o','markersize',5,'markerfacecolor','none','markeredgecolor','k');
            catch
                d = plot(X_Y_unq(subset,1),X_Y_unq(subset,2), '.','markersize',16,'markerfacecolor','none','markeredgecolor','k');
            end
            
        end
        %   plot([0 20],[0 20],'k--');
        hold off
        %axis equal
        set(gca,'fontsize',fontsize*.8);
        grid on
        cbar = colorbar;
        set(cbar,'fontsize',fontsize*.8);
        %xlim([0 10]);
        % ylim([0 1]);
        %  caxis([0.6 1]);
        % cbl = cblabel('swimming efficiency / best efficiency','fontsize',fontsize);
       
        % l = legend([d1 d2],'tail crudely optimized','tail not optimized');
        % set(l,'position',[0.58521      0.90343      0.23778     0.067114]);
        title(depvar,'interpreter','none')
        % cbl = cblabel('swimming speed (\mum/s)','fontsize',fontsize);
        % cbl = cblabel('motor frequency (rev/s)','fontsize',fontsize);
        %  cbl = cblabel('# tail wavelengths','fontsize',fontsize);
        %cbl = cblabel('tail amplitude (\mum)','fontsize',fontsize);
        %    cbl = cblabel('tail wavelength (\mum)','fontsize',fontsize);
        %  cbl = cblabel('# tails tested','fontsize',fontsize);
        %         best_line.AR1 = X(1,:);  %all AR1s
        %         best_line.AR2 = NaN(size(X(1,:)));  %best corresponding AR2s to be filled in below
        %         for i = 1:size(X,2)
        %             [~, max_eff_ind] = max(Z(:,i));
        %             best_line.AR2(i) = Y(max_eff_ind,1);
        %         end
        %         hold on
        %         pl = plot(best_line.AR1,best_line.AR2,'k-');
        %         hold off
        
        drawnow
%         saveas(gcf,[outfolder,prefixes{fi},'_',depvar,'_new.png'])
        
%         saveas(gcf,[outfolder,depvar,'_new.png'])
        
%            Depvars = {'Power_eff','Dm', 'temporal_SN','temporal_ability','fore_aft_SN','fore_aft_ability','tau_a','lambda_a','error_angle'};
%        
%            if ismember(dv,1:6)
%                Pareto_data.(depvar).Z = Z;
%                Pareto_data.(depvar).F = F;
%            end
        %%
    end  %depvars
    %%
    
end  %files

sdfsdfew
%% plot of % improvement of best curved vs straight for each AR1
% nearest neighbor interpolation above doesn't work for this - change to
% linear!


   method = 'linear';   
        [AR1_AR2,inds] = unique([ [Results.AR1]' [Results.AR2]' ],'rows');
        AR1_range = max(AR1_AR2(:,1)) - min(AR1_AR2(:,1));
        AR2_range = max(AR1_AR2(:,2)) - min(AR1_AR2(:,2));
        standardized_X_subset = (  X_Y_unq(subset,1) - min(AR1_AR2(:,1))  ) / AR1_range;
        standardized_Y_subset = (  X_Y_unq(subset,2) - min(AR1_AR2(:,2))  ) / AR2_range;
        F = scatteredInterpolant(standardized_X_subset, standardized_Y_subset, dep_var(subset),method,'none');
        standardized_X = (  X - min(AR1_AR2(:,1))  ) / AR1_range;
        standardized_Y = (  Y - min(AR1_AR2(:,2))  ) / AR2_range;
        Z = F(standardized_X,standardized_Y);


%plot best AR2 for each AR1 curve
xtemp = X(:); ytemp = Y(:); ztemp = Z(:);
unq_AR1s = unique(xtemp);  best_AR2s = NaN(size(unq_AR1s));  improvements = best_AR2s;
for i = 1:length(unq_AR1s)
    AR2_temp = ytemp( xtemp == unq_AR1s(i));  eff_temp = ztemp( xtemp == unq_AR1s(i));
    best_AR2s(i) = AR2_temp(eff_temp == max(eff_temp));
    improvements(i) = (max(eff_temp) - eff_temp( AR2_temp == 0) ) / eff_temp( AR2_temp == 0) * 100;
end

figure(581);  best_AR2_line = plot(unq_AR1s,best_AR2s,'ro--');  xlim([1 12]);  grid on;
xlabel('AR_1','fontsize',fontsize);  ylabel('best AR_2','fontsize',fontsize);
xlim([1 Inf])
set(gca,'fontsize',fontsize-4);
set(gca,'XTickLabel',[1 2:12]);
title('Best AR_2 for a given AR_1','fontsize',fontsize)



figure(582);  plot(unq_AR1s,improvements,'-','linewidth',2);  grid on
xlabel('AR_1','fontsize',fontsize);  ylabel('% improvement of curved vs straight','fontsize',fontsize);
xlim([1 Inf])
set(gca,'fontsize',fontsize-4);
set(gca,'XTickLabel',[1 2:12]);
title('Power eff improvement by being curved','fontsize',fontsize)


%%

varlist = {'AR1s','AR2s','amps','lambdas','nlambdas'};
% look at best choices of two vars, given each unique combo of all
% other vars
variable1 = 'AR1s';  variable2 = 'AR2s';
var1 = eff.curved_rod.(variable1);   var2 = eff.curved_rod.(variable2);
%var1 = eff.curved_rod.amps;  var2 = eff.curved_rod.lambdas;
others = setdiff(varlist,{variable1, variable2});

othervars = [];
for i = 1:length(others)
    othervars = [othervars  eff.curved_rod.(others{i})'];  %two vars to use for X and Y axes
end

[othervars_unq,inds,inds2] = unique(othervars,'rows');
% best.effs = NaN(size(othervars_unq,1),1);
fields = fieldnames(eff.curved_rod);
for f = 1:length(fields)
    field = fields{f};
    best.(field) = NaN(size(othervars_unq,1),1);
end
best.ind = NaN(size(othervars_unq,1),1);
best.npts = NaN(size(othervars_unq,1),1);
% best_var1s = best_depvars;
% best_var2s = best_depvars;

for i = 1:size(othervars_unq,1) %each unique combo of other vars
    othervar = othervars_unq(i,:);
    logical_temp = true(size(othervars,1),1);
    for j = 1:size(othervars,2)  %each var type, usually 3
        logical_temp = logical_temp & othervars(:,j) == othervar(j);  %incremental AND aggregation - if any is false, it goes false
    end
    inds_temp = find(logical_temp);
    effs = eff.curved_rod.power_effs(inds_temp);
    best_ind = inds_temp( effs == max(effs));  %index of best run within this subset
    % best_depvars(i) = eff.curved_rod.power_effs(best_ind);
    for f = 1:length(fields)
        best.(fields{f})(i) = eff.curved_rod.(fields{f})(best_ind);
    end
    best.ind(i) = best_ind;
    best.npts(i) = length(inds_temp);
    
end

%%
figure;  %plot(best.AR1s,best.AR2s,'o','markerfacecolor','k'); grid on
clear text
minpts = 5;
colorlo = 0.65;  colorhi = 1;
[colors,colorinds] = colordata(50,'parula',[colorlo colorhi],best.power_effs(best.npts > minpts) / max(best.power_effs));

x = best.AR1s(best.npts > minpts);  y = best.AR2s(best.npts > minpts);
for i = 1:length(x)
    plot(x(i),y(i),'o','markerfacecolor',colors(colorinds(i),:),'markeredgecolor',colors(colorinds(i),:),'markersize',6)
    if sum( x == x(i) & y == y(i)) > 1
        text(x(i)+0.01,y(i)+0.01,num2str(sum( x == x(i) & y == y(i))));
    end
    hold on
end
hold off

colorbar
caxis([colorlo colorhi]);
colormap(colors);
grid on
xlim([1 12]);  ylim([0 1]);
xlabel('AR1'); ylabel('AR2');
title('Best bodies for each unique tail tested')
cblabel('relative efficiency','fontsize',12)
%%
results.power_effs = NaN(length(AR1s),length(AR2s),length(amps(1,:)),length(lambdas(1,:)),length(nlambdas(1,:)));
results.speeds = results.power_effs;
results.freqs = results.speeds;
fields = fieldnames(results);

for i = 1:length(eff.curved_rod.AR1s)
    %     ind = find(AR1s == eff.curved_rod.AR1s(i),AR2s == eff.curved_rod.AR2s(i),amps == eff.curved_rod.amps(i),lambdas == eff.curved_rod.lambdas(i),nlambdas == eff.curved_rod.nlambdas(i));
    %     if length(ind) > 1
    %         stopafra
    %     end
    for f = 1:length(fields)
        results.(fields{f})(AR1s == eff.curved_rod.AR1s(i),AR2s == eff.curved_rod.AR2s(i),amps(1,:) == eff.curved_rod.amps(i),lambdas(1,:) == eff.curved_rod.lambdas(i),nlambdas(1,:) == eff.curved_rod.nlambdas(i)) = eff.curved_rod.(fields{f})(i);
    end
end

%%
min_tails = 29;
for i = 1:size(results.power_effs,1) %AR1
    for j = 1:size(results.power_effs,2) %AR2
        temp = ~isnan(results.power_effs(i,j,:,:,:));
        if sum(temp(:)) >= min_tails
            disp(['AR1 = ',num2str(AR1s(i)),'   AR2 = ',num2str(AR2s(j)),'   has ',num2str(sum(temp(:))),' tails']);
        end
    end
end
%%
colors = distinguishable_colors(20);
AR1 = 4.5;  AR2 = 0.625;
figure(45)
c = 0; clear legends
for i = 1:size(lambdas,2)
    for j = 1:size(nlambdas,2)
        
        if sum(~isnan(squeeze(results.power_effs(AR1s == AR1, AR2s == AR2, :,i,j)))) >= 0
            c = c+1;
            temp = squeeze(results.power_effs(AR1s == AR1, AR2s == AR2, :,i,j));
            x = amps(1,~isnan(temp));  y = temp(~isnan(temp));
            plot(x,y/max(results.power_effs(:)),'-o','markeredgecolor',colors(c,:),'markerfacecolor',colors(c,:),'color',colors(c,:));
            legends{c} = ['\lambda = ',num2str(lambdas(i)),',     # \lambda = ',num2str(nlambdas(j))];
            hold on
            grid on
            
        end
    end
end
hold off
xlabel('tail amplitude (\mum)','fontsize',fontsize);
ylabel('swimming efficiency / best efficiency','fontsize',fontsize);
legend(legends,'fontsize',fontsize-8);
title(['AR_1 = ',num2str(AR1),'       AR_2 = ',num2str(AR2)],'fontsize',fontsize);
set(gca,'fontsize',fontsize*.8);


%%
%%
sweeps_already_done = [eff.curved_rod.AR1s' eff.curved_rod.AR2s' eff.curved_rod.amps' eff.curved_rod.lambdas' eff.curved_rod.nlambdas'];
%%



AR1s = [1.75  1.75     2      3.5       5        6        11     12];
AR2s = [0.025 0.075   0.025   0.25     0.94      0.95     0.9    0.9];



method = 'natural';
method = 'nearest';
method = 'linear';

dep_var = best_amps;
F_amp = scatteredInterpolant(X_Y_unq(subset,1), X_Y_unq(subset,2), dep_var(subset),method,'none');
dep_var = best_lambdas;
F_lambda = scatteredInterpolant(X_Y_unq(subset,1), X_Y_unq(subset,2), dep_var(subset),method,'none');
dep_var = best_nlambdas;
F_nlambda = scatteredInterpolant(X_Y_unq(subset,1), X_Y_unq(subset,2), dep_var(subset),method,'none');

clear amps lambdas nlambdas
for i = 1:length(AR1s)
    amps(i) = F_amp(AR1s(i),AR2s(i));
    lambdas(i) = F_lambda(AR1s(i),AR2s(i));
    nlambdas(i) = F_nlambda(AR1s(i),AR2s(i));
end




%%

results_inds = NaN(size(X_Y_unq,1),1);
for i = 1:size(X_Y_unq,1)
    temp = find( [Results.AR1] == X_Y_unq(i,1) & [Results.AR2] == X_Y_unq(i,2) & [Results.amp] == best_amps(i) & [Results.lambda] == best_lambdas(i) & [Results.nlambda] == best_nlambdas(i));
    if length(temp) > 1
        temp
    end
    results_inds(i) = temp(end);
    % n = 2;
    % results_inds(i) = find(  roundn([Results.lambda],-n) == roundn(best_lambdas(i),-n) & roundn([Results.nlambda],-n) == roundn(best_nlambdas(i),-n)  );
end


%% straight rod line graph

straight_inds = find(X_Y_unq(:,2) == 0);
straight_bodies = X_Y_unq(straight_inds,:);
straight_dep_vars = dep_var(straight_inds);

figure(932)
plot(straight_bodies(:,1),straight_dep_vars,'o-','markerfacecolor','b');
grid on
xlabel('AR_1','fontsize',fontsize-3)
ylabel('Normalized Power eff','fontsize',fontsize-3)
title('Straight rod Power eff','fontsize',fontsize)
xlim([1 Inf]);
set(gca,'fontsize',fontsize-5)
hold on
plot([1.67 1.67],[0.55 1],'--r');
hold off
legend('this study','Shum et al ellipsoid optimum','location','best')

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

     AR1_new =      [1.5    2       2.5      3      3.5     4            5     5.5     6      6.5   7      7.5    8      8.5     9     9.5    10     10.5    11       11.5       12 ...
            7.5  8     8.5      9    9.5   10     10.5  11    11.5   12];
        AR2_new = [0.475  0.635   0.795    0.905  0.945   0.965        0.98  0.985   0.99   0.99  0.99   0.995  0.995  0.995   0.995 0.995  0.995  0.995   0.995    0.995      0.995 ...
            0.95   0.95  0.95  0.95  0.95    0.95  0.95  0.95  0.95  0.95];
       
        guess0 = [0.4 3 1.3];  %if we have no idea, use this guess
        
X_Y_unq_temp = [ X_Y_unq;  [AR1_new' AR2_new'] ];
best_depvars_temp = [ best_depvars; NaN(length(AR1_new),1)];
best_amps_temp = [ best_amps; NaN(length(AR1_new),1)];
best_lambdas_temp = [ best_lambdas; NaN(length(AR1_new),1)];
best_nlambdas_temp = [ best_nlambdas; NaN(length(AR1_new),1)];

% [~, AR1_inds] = sort(X_Y_unq(:,1));  bodies_AR1 = X_Y_unq(AR1_inds,:);  effs_AR1 = best_effs(AR1_inds);  amps_AR1 = best_amps(AR1_inds);  lambdas_AR1 = best_lambdas(AR1_inds);  nlambdas_AR1 = best_nlambdas(AR1_inds);
% [~, AR2_inds] = sort(X_Y_unq(:,2));  bodies_AR2 = X_Y_unq(AR2_inds,:);  effs_AR2 = best_effs(AR2_inds);  amps_AR2 = best_amps(AR2_inds);  lambdas_AR2 = best_lambdas(AR2_inds);  nlambdas_AR2 = best_nlambdas(AR2_inds);

[~, AR1_inds] = sort(X_Y_unq_temp(:,1));  bodies_AR1 = X_Y_unq_temp(AR1_inds,:);  effs_AR1 = best_depvars_temp(AR1_inds);  amps_AR1 = best_amps_temp(AR1_inds);  lambdas_AR1 = best_lambdas_temp(AR1_inds);  nlambdas_AR1 = best_nlambdas_temp(AR1_inds);
[~, AR2_inds] = sort(X_Y_unq_temp(:,2));  bodies_AR2 = X_Y_unq_temp(AR2_inds,:);  effs_AR2 = best_depvars_temp(AR2_inds);  amps_AR2 = best_amps_temp(AR2_inds);  lambdas_AR2 = best_lambdas_temp(AR2_inds);  nlambdas_AR2 = best_nlambdas_temp(AR2_inds);



clear max_der mean_der n_higher_effs n_neighbors
for i = 1:length(best_depvars_temp)
    body = X_Y_unq_temp(i,:);
    eff = best_depvars_temp(i);
    AR1_lo = find( bodies_AR1(:,1) < body(1) &  bodies_AR1(:,2) == body(2),1,'last');
    if isempty(AR1_lo)
        der_AR1_lo = NaN;
        eff_AR1_lo = NaN;
        amp_AR1_lo = NaN;  lambda_AR1_lo = NaN;  nlambda_AR1_lo = NaN;
    else
        der_AR1_lo = ( effs_AR1(bodies_AR1(:,1) == body(1) & bodies_AR1(:,2) == body(2))  -  effs_AR1(AR1_lo) )   / ( bodies_AR1(bodies_AR1(:,1) == body(1) & bodies_AR1(:,2) == body(2), 1)  -  bodies_AR1(AR1_lo , 1)  )    *   range(bodies_AR1(:,1))  /  effs_AR1(bodies_AR1(:,1) == body(1) & bodies_AR1(:,2) == body(2)) ;
        eff_AR1_lo = effs_AR1(AR1_lo);
        amp_AR1_lo = amps_AR1(AR1_lo);  lambda_AR1_lo = lambdas_AR1(AR1_lo);  nlambda_AR1_lo = nlambdas_AR1(AR1_lo);
    end
    
    AR1_hi = find( bodies_AR1(:,1) > body(1) &  bodies_AR1(:,2) == body(2),1,'first');
    if isempty(AR1_hi)
        der_AR1_hi = NaN;
        eff_AR1_hi = NaN;
        amp_AR1_hi = NaN;  lambda_AR1_hi = NaN;  nlambda_AR1_hi = NaN;
    else
        der_AR1_hi = (  effs_AR1(AR1_hi)  -  effs_AR1(bodies_AR1(:,1) == body(1) & bodies_AR1(:,2) == body(2)) )   / ( bodies_AR1(AR1_hi , 1)  -  bodies_AR1(bodies_AR1(:,1) == body(1) & bodies_AR1(:,2) == body(2), 1) )    *   range(bodies_AR1(:,1))  /  effs_AR1(bodies_AR1(:,1) == body(1) & bodies_AR1(:,2) == body(2)) ;
        eff_AR1_hi = effs_AR1(AR1_hi);
        amp_AR1_hi = amps_AR1(AR1_hi);  lambda_AR1_hi = lambdas_AR1(AR1_hi);  nlambda_AR1_hi = nlambdas_AR1(AR1_hi);
    end
    
    AR2_lo = find( bodies_AR2(:,2) < body(2) &  bodies_AR2(:,1) == body(1),1,'last');
    if isempty(AR2_lo)
        der_AR2_lo = NaN;
        eff_AR2_lo = NaN;
        amp_AR2_lo = NaN;  lambda_AR2_lo = NaN;  nlambda_AR2_lo = NaN;
    else
        der_AR2_lo = ( effs_AR2(bodies_AR2(:,2) == body(2) & bodies_AR2(:,1) == body(1))  -  effs_AR2(AR2_lo) )   / ( bodies_AR2( bodies_AR2(:,1) == body(1) & bodies_AR2(:,2) == body(2), 2)  -  bodies_AR2(AR2_lo , 2)  )    *   range(bodies_AR2(:,2))  /  effs_AR2(bodies_AR2(:,1) == body(1) & bodies_AR2(:,2) == body(2)) ;
        eff_AR2_lo = effs_AR2(AR2_lo);
        amp_AR2_lo = amps_AR2(AR2_lo);  lambda_AR2_lo = lambdas_AR2(AR2_lo);  nlambda_AR2_lo = nlambdas_AR2(AR2_lo);
    end
    
    AR2_hi = find( bodies_AR2(:,2) > body(2) &  bodies_AR2(:,1) == body(1),1,'first');
    if isempty(AR2_hi)
        der_AR2_hi = NaN;
        eff_AR2_hi = NaN;
        amp_AR2_hi = NaN;  lambda_AR2_hi = NaN;  nlambda_AR2_hi = NaN;
    else
        der_AR2_hi = ( effs_AR2(AR2_hi)  -  effs_AR2(bodies_AR2(:,2) == body(2) & bodies_AR2(:,1) == body(1)) )   / (bodies_AR2(AR2_hi , 2)   -  bodies_AR2( bodies_AR2(:,1) == body(1) & bodies_AR2(:,2) == body(2), 2) )    *   range(bodies_AR2(:,2))  /  effs_AR2(bodies_AR2(:,1) == body(1) & bodies_AR2(:,2) == body(2)) ;
        eff_AR2_hi =  effs_AR2(AR2_hi);
        amp_AR2_hi = amps_AR2(AR2_hi);  lambda_AR2_hi = lambdas_AR2(AR2_hi);  nlambda_AR2_hi = nlambdas_AR2(AR2_hi);
    end
    
    max_der(i) = max(abs([der_AR1_lo  der_AR1_hi  der_AR2_lo  der_AR2_hi]));
    mean_der(i) = nanmean(abs([der_AR1_lo  der_AR1_hi  der_AR2_lo  der_AR2_hi]));
    n_neighbors(i) = ~isnan(eff_AR1_lo) + ~isnan(eff_AR1_hi) + ~isnan(eff_AR2_lo) + ~isnan(eff_AR2_hi);
    n_higher_effs(i) = (eff_AR1_lo > eff) + (eff_AR1_hi > eff) + (eff_AR2_lo > eff) + (eff_AR2_hi > eff);
    
    amp_temp = [amp_AR1_lo amp_AR1_hi amp_AR2_lo amp_AR2_hi];
    lambda_temp = [lambda_AR1_lo lambda_AR1_hi lambda_AR2_lo lambda_AR2_hi];
    nlambda_temp = [nlambda_AR1_lo nlambda_AR1_hi nlambda_AR2_lo nlambda_AR2_hi];
    eff_temp = [eff_AR1_lo eff_AR1_hi eff_AR2_lo eff_AR2_hi];
    
    if isnan(max(eff_temp))
        guesses(i,:) = guess0;
    else
        guesses(i,:) = [amp_temp(eff_temp == max(eff_temp)) lambda_temp(eff_temp == max(eff_temp)) nlambda_temp(eff_temp == max(eff_temp))];
    end
    % shat = [eff       effs_AR1(bodies_AR1(:,1) == body(1) & bodies_AR1(:,2) == body(2))     effs_AR1(bodies_AR1(:,1) == body(1) & bodies_AR1(:,2) == body(2))      effs_AR2(bodies_AR2(:,2) == body(2) & bodies_AR2(:,1) == body(1))     effs_AR2(bodies_AR2(:,2) == body(2) & bodies_AR2(:,1) == body(1))     ]
    %     if length(unique(shat)) > 1
    %         stopa
    %     end
    
end

% [sorted_max_der, inds] = sort(max_der,'descend');


% mean_der( ~  (n_neighbors >= 3 & n_higher_effs >= 2) | (n_neighbors == 2 & n_higher_effs == 2)    ) = NaN;

% excludes = (X_Y_unq(:,1) >= 4 & X_Y_unq(:,1) <= 12 & X_Y_unq(:,2) <= 0.2 ) | X_Y_unq(:,2) >= 0.9 ;

% mean_der(excludes) = NaN;

% [sorted_mean_der, inds] = nansort(mean_der,'descend');

[sorted_mean_der, inds] = sort(mean_der,'descend');
% n_neighbors = n_neighbors(inds);  n_higher_effs = n_higher_effs(inds);
% X_Y_unq_sort = X_Y_unq(inds,:);



%inds = inds( (n_neighbors >= 3 & n_higher_effs >= 2) | (n_neighbors == 2 & n_higher_effs == 2)    );

bodies = X_Y_unq_temp(inds,:);
%guesses = [ best_amps(inds) best_lambdas(inds) best_nlambdas(inds) ];
guesses = guesses(inds,:);

% bodies =  [bodies;    4.75 0.625;  5.25 0.625; 4.75 0.65;  5.25 0.65;  4.75 0.675;  5.25 0.675;  ];
% guesses = [guesses;  [ repmat( best_amps(bodies(:,1)==5 & bodies(:,2)==0.65), 6, 1)  repmat( best_lambdas(bodies(:,1)==5 & bodies(:,2)==0.65), 6, 1)  repmat( best_nlambdas(bodies(:,1)==5 & bodies(:,2)==0.65), 6, 1)  ]  ];

bodies0 = bodies;  guesses0 = guesses;
%%
for n = 1:4
% n = 4;
bodies = bodies0(n:4:end,:);
guesses = guesses0(n:4:end,:);
save(['E:\Hull\sweep_CFD0',num2str(n),'.mat'],'bodies','guesses');
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

%%


AR1 = unique(X_Y_unq(:,1));  AR2 = unique(X_Y_unq(:,2));

for i = 1:length(AR1)
    
    inds = X_Y_unq(:,1) == AR1(i);
    AR2_plot = X_Y_unq(inds,2);
    depvar_plot = best_depvars(inds);
    
    [AR2_plot,inds] = sort(AR2_plot);
    depvar_plot = depvar_plot(inds);
    
    figure(84)
    plot(AR2_plot,depvar_plot,'o');  
    grid on
    title(['AR1 = ',num2str(AR1(i))]);
    
    pause
    
end