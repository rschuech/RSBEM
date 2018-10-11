function [s,e, varargout] = plot_mesh(Mesh,varargin)
%outputs surface handles (s) and edge handles (e)

% for plotting mesh - how much to refine in order to resolve curvature of
% quadratic triangles (since Matlab can only plot straight lines, flat
% patches)
% for i = 1:length(Mesh)
%     if Mesh(i).n_elem < 500
%         n_refines(i) = 3;
%     elseif Mesh(i).n_elem < 1000
%         n_refines(i) = 2;
%     elseif Mesh(i).n_elem < 1500
%         n_refines(i) = 2;
%     else
%         n_refines(i) = 2;
%     end
% end

if isempty(varargin)
    n_refines = 2;
else
    n_refines = varargin{1};
end

colors = {'g','r','r','b','c','c'};
facealpha = 1;
edgealpha = 0.9;
sidepts = 10;  %number of points along triangle edges was 7
% factors = [ 0.05 0.5 0.1 0.1];
factors = [ 0.2 0.75 0.15 0.15];
factors = [1 1 1 1 1 1 1];

if nargin >= 3
    Surface_vis = refine_vis_surface(Mesh,n_refines, varargin{2});
    for i = 1:length(Surface_vis)
        if length(varargin) <= 2
        color_limits(i,:) = [min(Surface_vis(i).surfun)  max(Surface_vis(i).surfun) * factors(i)];
        end
    end
    if length(varargin) <= 2
%     color_limits(3:4,:) = repmat( mean(color_limits(3:4,:),1) , 2, 1);
% color_limits = repmat( [min(color_limits(:)) max(color_limits(:))] , size(color_limits,1),1);
    else
        color_limits = varargin{3};
    end
else
    Surface_vis = refine_vis_surface(Mesh,n_refines);
end

Edge_vis = refine_vis_edges(Mesh,sidepts);


for i = 1:length(Surface_vis)
    if nargin >= 3
       [colors{i}] = colordata(500,'jet',color_limits(i,:),Surface_vis(i).surfun);
%          axh(i) = axes;
%          s(i) = patch(axh(i), 'faces',Surface_vis(i).elems,'vertices',Surface_vis(i).verts,'edgecolor','none','facecolor','interp','FaceVertexCData',Surface_vis(i).surfun,'facealpha',facealpha);
   
         s(i) = patch( 'faces',Surface_vis(i).elems,'vertices',Surface_vis(i).verts,'edgecolor','none','facecolor','interp','FaceVertexCData',colors{i},'facealpha',facealpha);
    
    else
%         axh = gca;
        s(i) = patch('faces',Surface_vis(i).elems,'vertices',Surface_vis(i).verts,'edgecolor','none','facecolor',colors{i},'facealpha',facealpha);
    end
    
    hold on
    
     if nargin < 3
    e(i) = patch('faces',Edge_vis(i).elems,'vertices',Edge_vis(i).verts,'facecolor','none','edgecolor','k','edgealpha',edgealpha);
     else
         e(i) = NaN;
     end
     
    axis equal; 
    
%     if i > 1 && nargin >= 3
%         axh(i).Visible = 'off';
%         axh(i).XTick = [];
%         axh(i).YTick = [];
%         axh(i).ZTick = [];
%     end
    
end
hold off
grid on;  box on;

h = zoom(gcf);
setAxes3DPanAndZoomStyle(h,gca,'camera');


% xlabel(axh(1), 'x');  ylabel(axh(1),'y');  zlabel(axh(1),'z');
xlabel( 'x');  ylabel('y');  zlabel('z');
% xtemp = [axh.XLim];  ytemp = [axh.YLim];  ztemp = [axh.ZLim];
% lims(:,:,1) = [xtemp(1:2:end)' xtemp(2:2:end)'; ] ;
% lims(:,:,2) = [ytemp(1:2:end)' ytemp(2:2:end)'; ] ;
% lims(:,:,3) = [ztemp(1:2:end)' ztemp(2:2:end)'; ] ;

% axlims = [squeeze(min(lims(:,1,:),[],1))    squeeze(max(lims(:,2,:),[],1)) ] ;
% 
% for i = 1:length(axh)
%     set(axh,'Xlim',axlims(1,:),'Ylim',axlims(2,:),'Zlim',axlims(3,:));
% end

% if nargin >= 3
% link = linkprop(axh, {'Box','CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle','Color','DataAspectRatio','FontSize','GridAlpha',...
%     'LineWidth','PlotBoxAspectRatio','Position','View','XLim','YLim','ZLim'});
% varargout = {link};
% end

%% Give each one its own colormap
% colormap(ax1,'hot')
% colormap(ax2,'cool')
%% Then add colorbars and get everything lined up
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);
% cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
% cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);

