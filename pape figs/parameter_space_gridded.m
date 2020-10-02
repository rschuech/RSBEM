
plot_type = 'all';  % 'all' for grid or 'selected' for a few separate figures or 'overlay' for just adding shapes to existing current fig
% plot_type = 'selected';
plot_type = 'overlay';
plot_subtype = 'pareto'; % 'pareto' for Pareto headliner, 'landscapes' for tasks:  from load_results, uses m = task index

switch plot_type
    case 'all'
        SF1 = 2:1:10;  SF2 = [0:0.1:0.8 0.9];  [X,Y] = ndgrid(SF1,SF2);
        bodies = unique([X(:) Y(:); 1 0; 1.125 0.1; 1.125 0.2; 1.125 0.3; 1.375 0.4; 1.625 0.5; 2.2 0.7;   3 0.8;     ],'rows');
        n_refines = 2;
    case 'selected'
        bodies = [1.4375 0;   10 0.15; 1 0; 4 0.7;  10 0;  6 0.8;  6 0.2;];  % should be 1.45 0 but I don't know where that mesh file is, prolly the servers
        % bodies = [5 0.7];
        n_refines = 2;
        
        %          SF1 = linspace(1,10,5);  SF2 = [0 0.25 0.5 0.75 0.95];  [X,Y] = ndgrid(SF1,SF2);
        %           bodies = unique([X(:) Y(:); 1 0;   1.05 0.25;       1.6 0.5;       2.4  0.75;        3.75 0.95; ],'rows');
        
        
        
    case 'overlay'
        
        switch plot_subtype
            case 'landscapes'
                
                SF1 = linspace(1,10,5);  SF2 = [0 0.25 0.5 0.75 0.95];  [X,Y] = ndgrid(SF1,SF2);
                bodies = unique([X(:) Y(:); ...
                    1 0;   1.3 0.25;   1.65 0.5;       2.45  0.75;        4 0.95; ...
                    ],'rows');   % for performance landscapes
                n_refines = 2;
                switch m
                    case 5
                        highlight_bodies = [1.46 0; 4 0.7]; %eff
                    case 6
                        highlight_bodies = [10 0.15]; %temporal
                    case {3, 22}
                        highlight_bodies = [1 0]; % tumbling, construction
                    case {1, 31, 32}
                        highlight_bodies = [10 0]; % uptake, tau_a body only, tau_a body+tail
                    case {33} 
                        highlight_bodies = [7 0.6]; % adj speed for SNR
                    otherwise
                        highlight_bodies = [];
                end
                
            case 'pareto'
                bodies = [];
                
%                 highlight_bodies = [      1        1.645        1.776        3.162        4.777        6.496        8.245       10.000        4.178;...
%                                           0        0.163        0.359        0.443        0.365        0.323        0.303        0.287        0.667]';
                 highlight_bodies = [      1         1.65         1.78        3.137        4.773        6.485        8.238       10.00     4.194;...
                                           0        0.166        0.361        0.443         0.37        0.322        0.298        0.292    0.665]';
                n_refines = 1;
        end
        
        is_highlighted = [ false(size(bodies,1),1);  true(size(highlight_bodies,1),1) ];
        
        bodies = [bodies; highlight_bodies];
          
          
  
end
fontsizes.axes = 20;   fontsizes.labels = 24;  fontsizes.title = 12;
% was 12 and 12 and 12
switch plot_type
    case {'all' , 'selected'}
        figure(827); clf; 
%         set(gcf,'Position',[ 672         237        1003         726]); 
          set(gcf,'Position',[   467         154        1022         730]);  % aspect ratio = 1.4, but mainly need to set for axis
    
        axes;
        set(gca,'Position',[ 0.13      0.13411        0.79      0.79]); % once fig is right aspect ratio, equal axes normalized width, height should make axes match fig aspect ratio
        set(gca,'PlotBoxAspectRatio',[   1      0.80206      0.80206]); % ALSO needed to reproduce aspect ratio....
    case 'overlay'
end





hold on;

switch plot_type
    case {'all' , 'selected'}
        grid on
        set(gca,'fontsize',fontsizes.axes)
        
        xlabel('$$ \mathcal{L} $$ (elongation)','interpreter','latex','fontsize',fontsizes.labels);
         ylabel('$$ \mathcal{K} $$ (curvature)','interpreter','latex','fontsize',fontsizes.labels);
        % title('Parameter Space','fontsize',fontsizes.title);
        
        xlim([1 10]);  ylim([0 1]);
        pbaspect([1      0.80206      0.80206]);
        
        
        set(gca,'YTick',[0:0.1:1]);
        set(gca,'XTick',[1:10]);
%         set(gca,'FontSize',13);
end

%%
if strcmp(plot_type,'overlay') && strcmp(plot_subtype,'pareto')
    ax = main_ax;
else
ax = gca;
end
meshfolder = 'C:\Hull\body meshes\';
clear axs bodplt surfaces edges
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
    
    if strcmp(plot_type,'overlay') && is_highlighted(b)
        ax_shape = axes('color','none','visible','off','clipping','off','position',[Xf Yf [0.019 0.019]*0.85  ]); % 0.7 t omake smaller
    else
        ax_shape = axes('color','none','visible','off','clipping','off','position',[Xf Yf [0.019 0.019]*0.85  ]); % 0.7 t omake smaller
    end
    %       set(ax_shape,'visible','on','box','on','xtick','','ytick','');
    axs(b) = ax_shape;
    
    
    orig_limits = [1 10; 0 1];
    shape_limits = [-20 20; -20 20] ;
    orig2shape.x = polyfit(orig_limits(1,:),shape_limits(1,:),1);
    orig2shape.y = polyfit(orig_limits(2,:),shape_limits(2,:),1);
    shape_coord = [polyval(orig2shape.x,bodies(b,1)) polyval(orig2shape.y,bodies(b,2)) 0]';
    %      Mesh = shiftMesh(Mesh,[shape_coord - Mesh(2).refpoints(:,1)]);
    
    if strcmp(plot_type,'selected')
        figure(500 + b); clf
    end
    [s,e] = plot_mesh(Mesh , n_refines);
    %     set(gca,'color','none','visible','off','clipping','off');
    switch plot_type
        case 'all'
            set(s,'FaceColor','g');
        case 'selected'
            set(s,'FaceColor',repmat(0.7,1,3));
        case 'overlay'
            if ~is_highlighted(b)
            set(s,'FaceColor',repmat(0,1,3),'facealpha',0.15);
            else
                switch plot_subtype
                    case 'landscapes'
                set(s,'FaceColor',repmat(0,1,3),'facealpha',1);  % black, for performance landscapes
                    case 'pareto'
                
                set(s,'FaceColor',repmat(0.,1,3),'facealpha',0.2);  % gray, for Pareto headliner
                end
            end
    end
    
    set(e,'edgealpha',0.0);
    
    surfaces(b) = s;
    edges(b) = e;
    %     xlim([Mesh(1).refpoints(1,1) Mesh(1).refpoints(1,1) + 1])
    %     ylim([Mesh(1).refpoints(2,1) Mesh(1).refpoints(2,1) + 1])
    
    switch plot_type
        case {'all' , 'overlay'}
            xlim([Mesh(1).Centroid(1) Mesh(1).Centroid(1) + 1])
            ylim([Mesh(1).Centroid(2) Mesh(1).Centroid(2) + 1])
        case 'selected'
            ylim([-1 3]);
            xlim([-4 5]);
            grid off
            set(gca,'color','none');
            axis off
            set(gcf,'Color','none');
    end
    
    zlim([-1 3.5])
    
    % hold on
    switch plot_type
        case {'all' , 'selected'}
            light
    end
    drawnow
    % pause
    if strcmp(plot_type,'selected')
        filename = ['C:\Hull\body shapes\',    'curved_rod_AR1_',num2str(bodies(b,1)),'_AR2_',num2str(bodies(b,2)),'.png'];
        %         print(filename,'-dpng');
        %         saveas(gcf,filename);
%         export_fig( filename, '-transparent');  % only export_fig seems able to make background transparent
    end
    
end

return

% this is the wrong function, should be export_fig.....
% exportfig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\gridded parameter space.png','format','png','color','rgb','resolution',800,'bounds','loose')
%%
for b = 1:length(surfaces)
    try, set(surfaces(b),'facecolor',repmat(0,1,3)); end
    
    try, set(surfaces(b),'facealpha',0.1); end
    
    
end