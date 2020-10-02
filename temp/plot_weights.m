


body = [7 0.5];  %test point to center weights function around
    

[AR1_AR2,inds] = unique([ [Results.AR1]' [Results.AR2]' ],'rows');
    AR1_range = max(AR1_AR2(:,1)) - min(AR1_AR2(:,1));
AR2_range = max(AR1_AR2(:,2)) - min(AR1_AR2(:,2));

       method = 'natural';
%         method = 'nearest';
        
         standardized_dists = sqrt(  1/AR1_range^2*(X - body(1) ).^2  +  1/AR2_range^2*(Y - body(2) ).^2  );

  %  ep = 2.5;
    rbf = @(r) 1./(1+(ep*r).^2); %see wikipedia on radial basis functions
    weights = rbf(standardized_dists);

        figure(787)

fontsize = 22;
        pcolor(X,Y,weights)
  colorbar
        shading interp
        xlabel('AR_1','fontsize',fontsize)
        ylabel('AR_2','fontsize',fontsize)

        hold on

 
           try
            d1 = plot(X_Y_unq(subset & isnan([next_iter.amp])',1),X_Y_unq(subset & isnan([next_iter.amp])',2), '.','markersize',8,'markerfacecolor','none','markeredgecolor','k');
            d2 = plot(X_Y_unq(subset & ~isnan([next_iter.amp])',1),X_Y_unq(subset & ~isnan([next_iter.amp])',2), 'o','markersize',5,'markerfacecolor','none','markeredgecolor','k');
           catch
               d = plot(X_Y_unq(subset,1),X_Y_unq(subset,2), '.','markersize',8,'markerfacecolor','none','markeredgecolor','k');
           end