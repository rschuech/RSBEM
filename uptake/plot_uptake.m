


  npts = 150;
[X,Y] = meshgrid(linspace(min(uptake.AR1),max(uptake.AR1),npts),linspace(min(uptake.AR2),max(uptake.AR2),npts));
    
uptake2 = uptake;
uptake2.flux(uptake2.AR1 == 1.8 & uptake2.AR2 == 0.45) = NaN;
uptake2.flux = uptake2.flux / uptake2.flux(uptake2.AR1 == 1 & uptake2.AR2 == 0);

AR1_AR2 = [ uptake2.AR1(~isnan(uptake2.flux))'  uptake2.AR2(~isnan(uptake2.flux))' ];
flux = uptake2.flux(~isnan(uptake2.flux))';

   method = 'linear';
    %               method = 'natural';
    
    AR1_range = max(AR1_AR2(:,1)) - min(AR1_AR2(:,1));
    AR2_range = max(AR1_AR2(:,2)) - min(AR1_AR2(:,2));
    limits = [min(AR1_AR2(:,1)) max(AR1_AR2(:,1)); min(AR1_AR2(:,2)) max(AR1_AR2(:,2))];
   
    [standardized] = standardize(AR1_AR2, limits);
    
    F = scatteredInterpolant(standardized(:,1), standardized(:,2), flux,method,'none');
    %         F = scatteredInterpolant(standardized_X_subset, standardized_Y_subset, dep_var(subset),method,'linear');
    
 
    [standardized_XY] = standardize([X(:) Y(:)], limits);
    standardized_X = NaN(size(X));  standardized_X(:) = standardized_XY(:,1);
    standardized_Y = NaN(size(Y));  standardized_Y(:) = standardized_XY(:,2);
    
    Z = F(standardized_X,standardized_Y);
   
    figure(354)
    pcolor(X,Y,Z);  shading interp;  hold on
      d = plot(AR1_AR2(:,1),AR1_AR2(:,2), '.','markersize',3,'markerfacecolor','none','markeredgecolor','k');
      hold off
      return
      %% incorporate uptake into Pareto data here
      % having just loaded uptake shat into variable [uptake]....
      Pareto_data.uptake.Z = uptake.Z;
      Pareto_data.uptake.F = uptake.F;
      
      for i = 1:size(AR1_AR2,1)
          ind =   find(ismember( [uptake.uptake2.AR1' uptake.uptake2.AR2']  , AR1_AR2(i,:), 'rows'));
          if numel(ind) == 1
              Pareto_data.uptake.depvar(i) = uptake.uptake2.flux(ind);
          else
              Pareto_data.uptake.depvar(i) = NaN;
          end
      end
      
      
      