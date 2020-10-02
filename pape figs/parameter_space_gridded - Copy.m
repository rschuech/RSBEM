

SF1 = 2:1:10;  SF2 = [0:0.1:0.8 0.9];  [X,Y] = ndgrid(SF1,SF2);
bodies = unique([X(:) Y(:); 1 0; 1.125 0.1; 1.125 0.2; 1.125 0.3; 1.375 0.4; 1.625 0.5; 2.2 0.7;   3 0.8;     ],'rows');

% bodies = [1.45 0;   10 0.15; 1 0; 4 0.7;];
% bodies = [5 0.7];

fontsizes.axes = 12;  fontsizes.labels = 12;  fontsizes.title = 12;

figure(827); clf; set(gcf,'Position',[672   543   678   420]); axes; hold on;


grid on
set(gca,'fontsize',fontsizes.axes)

xlabel('SF_1 (elongation)','fontsize',fontsizes.labels);
ylabel('SF_2 (curvature)','fontsize',fontsizes.labels);
% title('Parameter Space','fontsize',fontsizes.title);

xlim([1 10]);  ylim([0 1]);



set(gca,'YTick',[0:0.1:1]);
set(gca,'XTick',[1:10]);
set(gca,'FontSize',13);

ax = gca;
%%
meshfolder = 'C:\Hull\body meshes\';
clear axs bodplt
c = 0;
for b = 1:size(bodies,1)
    name = ['curved_rod_AR1_',num2str(bodies(b,1)),'_AR2_',num2str(bodies(b,2)),'.dat'];
    [Mesh, Metadata] = load_mesh([meshfolder,name],[],[]);
    if isempty(Mesh)
        continue
    end
    
    c = c + 1;
    
    bodplt(c,:) = bodies(b,:);
    
    [Mesh] = global_inds(Mesh);
    [Mesh] = renumber_Mesh(Mesh);
    
    clear temp_input
    temp_input.performance.nthreads = 8;
    temp_input.accuracy.integration_tol.area.abstol = 1E-6;
    temp_input.accuracy.integration_tol.area.reltol = 1000;
    temp_input.accuracy.integration_tol.area.maxevals = Inf;
    
    temp_input.accuracy.integration_tol.centroid.abstol = 1E-6;
    temp_input.accuracy.integration_tol.centroid.reltol = 1000;
    temp_input.accuracy.integration_tol.centroid.maxevals = Inf;
    
    temp_input.accuracy.integration_tol.volume.abstol = 1E-6;
    temp_input.accuracy.integration_tol.volume.reltol = 1000;
    temp_input.accuracy.integration_tol.volume.maxevals = Inf;
    
    [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
    
    if bodies(b,1) == 1 %sphere degenerate case
        Mesh(1).refpoints = [0 0 0]';   %center
    elseif bodies(b,2) == 0  %not a sphere, but a straight rod degenerate case
        Mesh(1).refpoints = [0 0 0]';   %center
    else  %actually a general curved rod
        Mesh(1).refpoints =  [2*Metadata(1).geom.radius_curv*sin(Metadata(1).geom.nturns / 2 *pi) * cos(Metadata(1).geom.nturns / 2 *pi); ...
            2*Metadata(1).geom.radius_curv*(sin(Metadata(1).geom.nturns / 2 *pi))^2;
            0 ];   %center of centerline
    end
    
    
    [Xf, Yf] = ds2nfu(ax,bodies(b,1) + 0., bodies(b,2) - 0.) ;
    ax_shape = axes('color','none','visible','off','clipping','off','position',[Xf Yf 0.019 0.019]);
    %       set(ax_shape,'visible','on','box','on','xtick','','ytick','');
    axs(b) = ax_shape;
    
    
    orig_limits = [1 10; 0 1];
    shape_limits = [-20 20; -20 20];
    orig2shape.x = polyfit(orig_limits(1,:),shape_limits(1,:),1);
    orig2shape.y = polyfit(orig_limits(2,:),shape_limits(2,:),1);
    shape_coord = [polyval(orig2shape.x,bodies(b,1)) polyval(orig2shape.y,bodies(b,2)) 0]';
    %      Mesh = shiftMesh(Mesh,[shape_coord - Mesh(2).refpoints(:,1)]);
    
%             figure;  
    [s,e] = plot_mesh(Mesh);
%     set(gca,'color','none','visible','off','clipping','off');
    
set(s,'FaceColor','g');

    set(e,'edgealpha',0.0);
    
    surfaces(b) = s;
    edges(b) = e;
%     xlim([Mesh(1).refpoints(1,1) Mesh(1).refpoints(1,1) + 1])
%     ylim([Mesh(1).refpoints(2,1) Mesh(1).refpoints(2,1) + 1])
        xlim([Mesh(1).Centroid(1) Mesh(1).Centroid(1) + 1])
    ylim([Mesh(1).Centroid(2) Mesh(1).Centroid(2) + 1])
    
    
    zlim([-1 3.5])
    
    % hold on
%     light
    drawnow
    % pause
end
%%
for b = 1:length(surfaces)
    try, set(surfaces(b),'facecolor',repmat(0,1,3)); end
    
     try, set(surfaces(b),'facealpha',0.1); end
     
     
end