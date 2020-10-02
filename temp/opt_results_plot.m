

dep = [Results.Power_eff];
amp = [Results.amp];  lambda = [Results.lambda];  nlambda = [Results.nlambda];

[dep,inds] = sort(dep,'descend');  amp = amp(inds);  lambda = lambda(inds);  nlambda = nlambda(inds);


if ~isempty(dep)
[colors,colorinds] = colordata(200,'jet',[min(dep) max(dep)],dep);
else
    colors = [];  colorinds = [];
end

markersize = 8;
if ~isempty(dep)
figure(830)
% subplot(1,3,1)
for i = 2:length(amp)
plot3(amp(i),lambda(i),nlambda(i),'o','markersize',markersize,'markerfacecolor',colors(colorinds(i),:),'markeredgecolor',colors(colorinds(i),:))
text(amp(i)+( max(amp)-min(amp))*0.001,lambda(i)+( max(lambda)-min(lambda))*0.001,nlambda(i)+( max(nlambda)-min(nlambda))*0.001,{num2str(dep(i)/max(dep))} );
hold on
end
plot3(amp(1),lambda(1),nlambda(1),'p','markersize',markersize+5,'markerfacecolor',colors(colorinds(1),:),'markeredgecolor',colors(colorinds(1),:))
%text(amp(1)+( max(amp)-min(amp))*0.02,lambda(1)+( max(lambda)-min(lambda))*0.02,nlambda(1),{num2str(dep(1)/max(dep)),num2str(reldiffs(1))} );
% plot3(pm_guess(1),pm_guess(2),pm_guess(3),'*','markersize',markersize+5,'markerfacecolor','k','markeredgecolor','k')
hold off
xlabel('amp'); ylabel('lambda'); zlabel('nlambda');
view([0 90]);
caxis([min(dep) max(dep)])
colormap(colors)
% title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))]});
end








.4486 - .4642

ans =

      -0.0156

3.24 - 3.16

ans =

         0.08

1.323 - 1.292

ans =

        0.031