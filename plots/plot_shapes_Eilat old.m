
dumpfolder = 'E:\Hull\select_dumps\';

dumps = dir([dumpfolder,'*.mat']);
dumps = {dumps.name};

for d = 1:length(dumps)-3
    d
    load([dumpfolder,dumps{d}]);
    Mesh.Volume
    figure(d+200)
    clf
    
    [s,e] = plot_mesh(Mesh);  light;
    set(e,'edgealpha',0.2);
       xlim([-9 3.25])
    ylim([-1 3.5])
    zlim([-1 3.5])
    
    drawnow
%     pause
    
end

%%

for d = 1:length(dumps)
   
   
    figure(d)
%     set(gcf,'renderer','painters')
%     xlim([-9 3.25])
%     ylim([-1 3.5])
    axis off
%      axis on
    set(gcf,'color','w')
   name = [ dumps{d}(1:end-4)  ];
      print('-dpng','-r500',[dumpfolder,name,'.png']);
%         print([dumpfolder,name,'.eps'],'-depsc','-tiff','-painters');
    
end

%%
fontsizes.labels = 32;
fontsizes.axes = 24;
fontsizes.title = 33;

bodies = [1 0; 10 0; 10 0.95; 1.125 0.35; 2 0.635; 3.5 0.945; 5.5 0.985; 5.5 0; 10 0.65; 6 0.65; 3.5 0.25; 5.5 0.15;   10 0.3;  6 0.35; 5.5 0.5; 3 0.5; 10 0.5;];


figure(592)

% box on
xlim([1 10]);  ylim([0 1]);
plot(bodies(:,1),bodies(:,2),'ko','markerfacecolor','k');
grid on
set(gca,'fontsize',fontsizes.axes)

xlabel('AR_1 (elongation)','fontsize',fontsizes.labels);
ylabel('AR_2 (curvature)','fontsize',fontsizes.labels);
title('Parameter Space','fontsize',fontsizes.title);



set(gcf,'Position',[    519          38        1132         971]);
% set(gcf,'PaperPositionMode','manual');
 set(gcf,'PaperPosition',[    -1.6458      0.44271       14.15       12.137  ]);
 set(gca,'XTick',[1:10]);
   print('-dpng','-r600',[dumpfolder,'parameter_space','.png']);
