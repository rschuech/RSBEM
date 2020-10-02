
AR1 = linspace(1,10,100);
V = 1;   sph_rad = (V*3/4/pi)^(1/3);
mu = 1E-9;

%% eff based on friction coefficients of ellipsoids ala Dusenbery
% b = c
% AR1 = a/b
% V = 4/3 pi a b^2
% 
% b = sqrt(3/4 V / pi / a) 
% 
% AR1 = a / [ sqrt(3/4 V / pi) 1/sqrt(a)  ]  =  a^(1.5) / sqrt(3/4 V / pi) 

a = (  AR1 .* sqrt(3/4 .* V ./ pi)  ).^(2/3);
b = a ./ AR1;
c = b;

clear f_translation f_rotation freeswim text
for i = 1:length(AR1)
    clear geom
    geom.a = a(i);  geom.b = b(i);  geom.c = c(i);
    [f_translation.ellipsoid(i,:), f_rotation.ellipsoid(i,:)] = ellipsoid_theoretical_friction_coeffs(geom, mu);
end
    

% P = f * U^2 so the U cancels out in eff
eff_ellipsoid = 6*pi*sph_rad*mu ./ f_translation.ellipsoid(:,1);

pp = spline(AR1,eff_ellipsoid);
[best_straight.AR1.ellipsoid,temp] = fminbnd(@(ar1) -ppval(pp,ar1),1,3);  best_straight.eff.ellipsoid = -temp;
%% load straight rod body only fcoeffs from simulations
capsule_data = load('C:\Hull\CFD03\dino_code\pape figs\straight_rod_body_only_friction_coeffs.mat');
f_translation.capsule = vertcat(capsule_data.F_coeffs.translation);
eff_capsule = 6*pi*sph_rad*mu ./ f_translation.capsule(:,1);  eff_capsule = eff_capsule / eff_capsule(1);  % normalize to sphere
pp = spline(capsule_data.AR1,eff_capsule);
[best_straight.AR1.capsule,temp] = fminbnd(@(ar1) -ppval(pp,ar1),1,3); best_straight.eff.capsule = -temp;
%% straight rod free swimming efficiencies
inds = find(Best.eff.body(:,2) == 0);
% ind = find(Best.eff.body(:,2) == 0 & Best.eff.body(:,1) == 1.45);
% inds = setdiff(inds,ind);
freeswim.eff.AR1 = Best.eff.body(inds,1);
freeswim.eff.metric = Best.eff.metric(inds);
[freeswim.eff.AR1,inds] = sort(freeswim.eff.AR1);
freeswim.eff.metric = freeswim.eff.metric(inds);  sph_eff = freeswim.eff.metric(1);  freeswim.eff.metric = freeswim.eff.metric / sph_eff;
pp = spline(freeswim.eff.AR1,freeswim.eff.metric);
[best_straight.AR1.freeswim,temp] = fminbnd(@(ar1) -ppval(pp,ar1),1,3);  best_straight.eff.freeswim = -temp;
%% make graph
best_straight.AR1.shum = 1.67;  best_straight.eff.shum = 6.43E-3 / sph_eff; % can't find sphere results in Shum et al, will trust that my sphere is same as his

fontsize = 16;  fontsize_annotations = 12;
figure(821); clf;  set(gcf,'Position',[ 479          96        1346         882]); set(gcf,'Color','w');
plot(AR1, eff_ellipsoid, 'b-','linewidth',1.5);
hold on
plot(capsule_data.AR1', eff_capsule, 'r-','linewidth',1.5,'Marker','.','MarkerSize',10);
plot(1.67,best_straight.eff.shum,'o','markerfacecolor',repmat(0.5,1,3),'MarkerEdgeColor',repmat(0.5,1,3),'markersize',7)
plot(freeswim.eff.AR1,freeswim.eff.metric,'k-','linewidth',1.5,'Marker','.','MarkerSize',10);

plot(repmat(best_straight.AR1.ellipsoid,1,2),[0.9 best_straight.eff.ellipsoid],'b:','linewidth',1.);
plot(repmat(best_straight.AR1.capsule,1,2),[0.875 best_straight.eff.capsule],'r:','linewidth',1.);

plot(repmat(best_straight.AR1.freeswim,1,2),[0.825 best_straight.eff.freeswim],'k:','linewidth',1.);


plot(repmat(best_straight.AR1.shum,1,2),[0.85 best_straight.eff.shum],':','linewidth',1.,'Color',repmat(0.5,1,3));
text(best_straight.AR1.ellipsoid,0.9-8E-3,['$$ \mathcal{L} $$ = ',num2str(best_straight.AR1.ellipsoid,3)],'interpreter','latex','color','b','fontsize',fontsize_annotations);
text(best_straight.AR1.capsule,0.875-8E-3,['$$ \mathcal{L} $$ = ',num2str(best_straight.AR1.capsule,3)],'interpreter','latex','color','r','fontsize',fontsize_annotations);
text(best_straight.AR1.shum,0.85-8E-3,['$$ \mathcal{L} $$ = ',num2str(best_straight.AR1.shum,3)],'interpreter','latex','color',repmat(0.5,1,3),'fontsize',fontsize_annotations);
text(best_straight.AR1.freeswim,0.825-8E-3,['$$ \mathcal{L} $$ = ',num2str(best_straight.AR1.freeswim,3)],'interpreter','latex','color','k','fontsize',fontsize_annotations);

% xtemp = linspace(1,10,5000);
% ytemp = fnval(pp,xtemp);
% plot(xtemp,ytemp,'b--');

hold off
grid off
legend('ellipsoid, body-only (i.e., Dusenbery)','capsule, body-only (this study)','ellipsoid, flagellum-driven (Shum et al)','capsule, flagellum-driven (this study)','Location','best');
xlabel('$$ \mathcal{L} $$ (elongation)','interpreter','latex');  ylabel('relative swimming efficiency');
set(gca,'FontSize',fontsize);
% return

% export_fig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\straight rod efficiency comparisons.png','-r600','-nocrop')


%% 


tau_a_rel = f_rotation.ellipsoid(:,2) ./ f_rotation.ellipsoid(1,2);  % tau_a normalized to tau_a of sphere is simply the ratio of rotational friction coeffs for either axis b or c (same for ellipsoids of rev), Dusenbery p. 111
f_a_rel = f_translation.ellipsoid(:,1) ./ f_translation.ellipsoid(1,1); % translational friction coeff along long axis normalized to that of a sphere
v_rel = f_a_rel.^(-1/2); % swimming speed relative to that of sphere (since v^2 is proportional to 1/f , Dusenbery p. 177)

Dm_rel = tau_a_rel ./ f_a_rel;  % motile diffusivity relative to sphere (Dusenbery p. 177)
temporal_SN_rel = tau_a_rel.^(3/2) ./ f_a_rel.^(1/2); % temporal S/N relative to sphere (Dusenbery p. 255)
fore_aft_SN_rel = (a'/sph_rad) .* tau_a_rel.^(1/2); % fore-aft S/N relative to sphere (Dusenbery p. 255)

inds = find(Best.Dm.body(:,2) == 0);
freeswim.Dm.AR1 = Best.Dm.body(inds,1);
freeswim.Dm.metric = Best.Dm.metric(inds);
[freeswim.Dm.AR1,inds] = sort(freeswim.Dm.AR1);
freeswim.Dm.metric = freeswim.Dm.metric(inds);  sph_Dm = freeswim.Dm.metric(1);  freeswim.Dm.metric = freeswim.Dm.metric / sph_Dm;

inds = find(Best.temporal.body(:,2) == 0);
freeswim.temporal.AR1 = Best.temporal.body(inds,1);
freeswim.temporal.metric = Best.temporal.metric(inds);
freeswim.temporal.v = Best.temporal.Adj_Speed(inds);
freeswim.temporal.tau_a = Best.temporal.tau_a(inds);
[freeswim.temporal.AR1,inds] = sort(freeswim.temporal.AR1);
freeswim.temporal.metric = freeswim.temporal.metric(inds);  sph_temporal = freeswim.temporal.metric(1);  freeswim.temporal.metric = freeswim.temporal.metric / sph_temporal;
freeswim.temporal.v = freeswim.temporal.v(inds);  sph_v = freeswim.temporal.v(1);  freeswim.temporal.v = freeswim.temporal.v / sph_v;
freeswim.temporal.tau_a = freeswim.temporal.tau_a(inds);  sph_tau_a = freeswim.temporal.tau_a(1);  freeswim.temporal.tau_a = freeswim.temporal.tau_a / sph_tau_a;

inds = find(Best.fore_aft.body(:,2) == 0);
freeswim.fore_aft.AR1 = Best.fore_aft.body(inds,1);
freeswim.fore_aft.metric = Best.fore_aft.metric(inds);
[freeswim.fore_aft.AR1,inds] = sort(freeswim.fore_aft.AR1);
freeswim.fore_aft.metric = freeswim.fore_aft.metric(inds);  sph_fore_aft = freeswim.fore_aft.metric(1);  freeswim.fore_aft.metric = freeswim.fore_aft.metric / sph_fore_aft;




figure(825);  set(gcf,'color','w');  clf
% h1 = plot(AR1,Dm_rel,'b:',freeswim.Dm.AR1,freeswim.Dm.metric,'r:');  hold on
% h2 = plot(AR1,fore_aft_SN_rel,'b--',freeswim.fore_aft.AR1,freeswim.fore_aft.metric,'r--');
h3 = plot(AR1,temporal_SN_rel,'k--',freeswim.temporal.AR1,freeswim.temporal.metric,'k-');  hold off;   set(h3,'linewidth',2)
% h = [h1; h2; h3];  set(h,'linewidth',2)


% plot(AR1,v_rel.^2,'b-',freeswim.temporal.AR1,freeswim.temporal.v.^2,'r-');   % Adj_Speed going with optimized temporal S/N is noisy!  but pretty small difference between Dusenbery and simulated relative speeds

% plot(AR1,tau_a_rel.^(3/2),'b-',freeswim.temporal.AR1,freeswim.temporal.tau_a.^(3/2),'r-'); % difference in SN is basically due to tau

grid on
ylim([1 45]);  xlim([1 Inf]);
xlabel('$$ \mathcal{L} $$ (elongation)','interpreter','latex');  ylabel('performance relative to spherical');
% legend('motile dispersal (Dusenbery)','motile dispersal (this study)','fore-aft S/N (Dusenbery)','fore-aft S/N (this study)','temporal S/N (Dusenbery)','temporal S/N (this study)','location','best');
legend('Dusenbery','this study','location','best');
set(gca,'FontSize',12);
set(gca,'YTick',[1 5:5:45]);

% export_fig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\straight rod SNR Dusenbery comparison.png','-r600','-nocrop')