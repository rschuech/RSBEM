


titles = {'Swimming Efficiency','Dispersal','Temporal SN','Temporal Chemotaxis','Fore-Aft SN','Fore-Aft Chemotaxis','Construction Ease'};
Depvars = {'Power_eff','Dm', 'temporal_SN','temporal_ability','fore_aft_SN','fore_aft_ability','construction'}; %,'tau_a','lambda_a','error_angle'};
Depvars = {'Power_eff','Dm', 'temporal_SN','fore_aft_SN','construction','tumbling','buckling_force_mean','buckling_force_max','buckling_torque_mean','buckling_torque_max'}; %,'tau_a','lambda_a','error_angle'};
titles = {'Swimming Efficiency','Dispersal','Temporal S/N','Fore-Aft S/N','Construction Ease','Tumbling Ease','Flagellar Buckling Force (mean)','Flagellar Buckling Force (max)','Flagellar Buckling Torque (mean)','Flagellar Buckling Torque (max)'};


titles = {'optimal flagellum amplitude','optimal flagellum wavelength','optimal flagellum # wavelengths'};
Fnames = {'amp','lambda','nlambda'};

dv = 3;  % temporal SN
% dv = 1;  % efficiency



depvar = Depvars{dv};

fontsize = 22;



npts = 1400;
npts = 150;
plot_datapts = true;



if ~strcmp(depvar,'construction')
    
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
    
else
    best_depvars =  Pareto_data.Construction_Ease.F(standardize(X_Y_unq, limits));
end

dep_vars = [best_amps best_lambdas best_nlambdas];

for ddv = 1:size(dep_vars,2)
    
    figure(ddv+900)
    
    
    best_depvars = dep_vars(:,ddv);
    
    [X,Y] = meshgrid(linspace(min(var1),max(var1),npts),linspace(min(var2),max(var2),npts));
    
    
    %         dep_var = best_depvars / max([Results.(depvar)]);  %normalize to global best performance
    
    [~,temp] = ismember(X_Y_unq,[1 0],'rows');  ind = find(temp);  % ind for sphere, probably first one...
    
    sphere_depvar = best_depvars(ind);
    
%     dep_var = best_depvars / sphere_depvar;  %normalize to sphere performance
    
    dep_var = best_depvars;
    
    subset = n_pts >= 1;
    %
    %         method = 'nearest';
    method = 'linear';
    %               method = 'natural';
    
    [AR1_AR2,inds] = unique([ [Results.AR1]' [Results.AR2]' ],'rows');
    
%     AR1_range = max(AR1_AR2(:,1)) - min(AR1_AR2(:,1));
%     AR2_range = max(AR1_AR2(:,2)) - min(AR1_AR2(:,2));
%     limits = [min(AR1_AR2(:,1)) max(AR1_AR2(:,1)); min(AR1_AR2(:,2)) max(AR1_AR2(:,2))];
     limits = [1 10; 0 1];
    %         standardized_X_subset = (  X_Y_unq(subset,1) - min(AR1_AR2(:,1))  ) / AR1_range;
    %         standardized_Y_subset = (  X_Y_unq(subset,2) - min(AR1_AR2(:,2))  ) / AR2_range;
    [standardized_subset] = standardize(X_Y_unq(subset,:), limits);
    
    F = scatteredInterpolant(standardized_subset(:,1), standardized_subset(:,2), dep_var(subset),method,'none');
    %         F = scatteredInterpolant(standardized_X_subset, standardized_Y_subset, dep_var(subset),method,'linear');
    F_tails.(Fnames{ddv}) = F;
    
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
    
    if ~strcmp(depvar,'construction')
        Pareto_data.(Depvars{dv}).Z = Z;
        Pareto_data.(Depvars{dv}).F = F;
        Pareto_data.(Depvars{dv}).depvar = best_depvars;
    end
    %
    
    fontsizes.labels = 30;
    fontsizes.axes = 22;
    fontsizes.title = 31;
    fontsizes.caxis = 26;
    
    %the switcheroo is used here so that standardized locations were used to do
    %interpolation, but orig locations are used for plot
    pcolor(X,Y,Z)
    %shading interp
    shading flat
    
    set(gca,'fontsize',fontsizes.axes);
    
    
    xlabel('AR_1 (elongation)','fontsize',fontsizes.labels)
    ylabel('AR_2 (curvature)','fontsize',fontsizes.labels)
    %
    %             xlabel('tail amplitude (\mum)','fontsize',fontsize)
    %         ylabel('tail wavelength (\mum)','fontsize',fontsize)
    %           ylabel('# tail wavelengths','fontsize',fontsize)
    %zlabel('power eff')
    hold on
    
    % plot single best pt
    xtemp = X_Y_unq(subset,1); ytemp = X_Y_unq(subset,2); ztemp =  dep_var(subset);
    %         best_pt =   plot(xtemp(ztemp==max(ztemp(:))),ytemp(ztemp==max(ztemp(:))),'kh','markerfacecolor','k','markersize',12);
    
    
    
    
    
    if plot_datapts
        
        d = plot(X_Y_unq(subset,1),X_Y_unq(subset,2), '.','markersize',10,'markerfacecolor','none','markeredgecolor','k');
%         switch Depvars{dv}
%             case 'Power_eff'
%                 starh = plot(6,0.65,'kp','markersize',20,'markerfacecolor','k');
%             case 'construction'
%                 starh = plot(1,0,'kp','markersize',20,'markerfacecolor','k');
%         end
        
    end
    
    hold off
    
    cbar = colorbar;
    
    
    % cbl = cblabel('swimming efficiency / best efficiency','fontsize',fontsize);
    
    % l = legend([d1 d2],'tail crudely optimized','tail not optimized');
    % set(l,'position',[0.58521      0.90343      0.23778     0.067114]);
    %         title(depvar,'interpreter','none')
    title(titles{ddv},'fontsize',fontsizes.title,'interpreter','none');
    % cbl = cblabel('swimming speed (\mum/s)','fontsize',fontsize);
    % cbl = cblabel('motor frequency (rev/s)','fontsize',fontsize);
    %  cbl = cblabel('# tail wavelengths','fontsize',fontsize);
    %cbl = cblabel('tail amplitude (\mum)','fontsize',fontsize);
    %    cbl = cblabel('tail wavelength (\mum)','fontsize',fontsize);
    
    
    
    
    
    hold on
    %         pl = plot(best_line.AR1,best_line.AR2,'k-','linewidth',3);
    
    
    clims = [ 0.8  1.6;       0.9  1.5;    0.9  1.1;          0.0  6.4;          0.6  1;             0.05  1;    0.9 1.4;                     0.9 1.4;                                                       0.95 12;             0.95 12;];
    cvals = {[0.85:0.01:0.9 1  1.2:0.05:1.6 1.475],[0.85:0.01:0.9 1  1.2:0.05:1.6 1.475],[1.05 1.2 1.5 1.8 2],[0.1 0.4 0.8 1 1.5 2 3 4 5 6],[0.65 0.7 0.75 0.8 0.9 0.95],[0.1:0.1:1],[ 0.95 0.975 1 1.125 1.2 1.25],     [0.95 1 1.1 1.2 1.3 1.4 1.5],    [1.1 1.5 3 2:2:10],    [1.1 1.5 3 2:2:10] };
    %caxis(clims(ddv,:));
    %   set(cbar,'fontsize',fontsizes.caxis);  % gets overidden by
    %   cblabel, apparently no way to fix
%     cbl = cblabel('performance relative to spherical body','fontsize',fontsizes.labels);
    
    hold on

%     if use_cvals
%     [C,ch] = contour(unique(X),unique(Y),Z,cvals{dv},'k--','linewidth',1);
%     else
        [C,ch] = contour(unique(X),unique(Y),Z,30,'k--','linewidth',1);
%     end
    
    
    %[C,ch] = contour(unique(X),unique(Y),Z,cvals{ddv},'k--','linewidth',1);
    %clabel(C,ch,'fontsize',fontsizes.axes,'LabelSpacing',1000);
    hold off
    % set(ch,'linewidth',2);
    set(gcf,'Position',[ 680          82        1051         896]);
    ylim([0 1]);
    set(gca,'XTick',[1:10])
    
    
    
    
    
    hold on
        % plot measurements
    try, delete(hmeas); end;
    for i = 1:length(plot_data)
        %      hmeas(i) = plot((plot_data(i).AR1),(plot_data(i).AR2),'bo','markerfacecolor','b','markersize',1);
      %  hmeas(i) = plot((plot_data(i).mean_unweighted(1)),(plot_data(i).mean_unweighted(2)),'ko','markerfacecolor','k','markersize',8);
    end
    hold off
    
    xlim([1 10]);
    
    drawnow
    %           print('-dpng','-r600',['E:\Hull\select_dumps\',titles{dv},'.png']);
    %         pause
    %         saveas(gcf,[outfolder,prefixes{fi},'_',depvar,'_new.png'])
    
    %         saveas(gcf,[outfolder,depvar,'_new.png'])
%  print('-dpng','-r300',['E:\Hull\graphs2\',titles{ddv},'.png']);
    
end  %ddv

