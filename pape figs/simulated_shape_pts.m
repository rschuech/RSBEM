figure(8293); set(gcf,'Color','w');
clf

plot(AR1_AR2_data.eff(:,1),AR1_AR2_data.eff(:,2), '.','markersize',4,'markerfacecolor','none','markeredgecolor',repmat(0,1,3));
hold on

plot(AR1_AR2_data.Dm(:,1),AR1_AR2_data.Dm(:,2), '.','markersize',11,'markerfacecolor','k','markeredgecolor','k');
% only using Dm not temporal since I've been marking some temporal points
% as bads and not loading them lately

% plot(AR1_AR2_data.uptake(:,1),AR1_AR2_data.uptake(:,2), '.','markersize',11,'markerfacecolor','k','markeredgecolor','k');


xlim([0.925 10.1]); ylim([-0.015 1]);

set(gca,'FontSize',16);

xlabel(gca,'$$ \mathcal{L} $$ (elongation)','interpreter','latex','fontsize',18);

ylabel(gca,'$$ \mathcal{K} $$ (elongation)','interpreter','latex','fontsize',18);

set(gcf,'Position',[   467         154        1022         730]);  % aspect ratio = 1.4, but mainly need to set for axis

set(gca,'position',[  0.132863531353356    0.174766424535667  0.79   0.79]);
set(gca,'PlotBoxAspectRatio',[   1      0.80206      0.80206]); % ALSO needed to reproduce aspect ratio....
