
plot_links = true;
color_links = true;

t = 0;
tfinal = sol.x(end);
%  tfinal = 0.14;
% tfinal = 0.54;
% tfinal = 2;


% Repulsion = calc_repulsive_forces(Mesh0, Network0, index_mapping, assembly_input.repulsion);

f_s = Network0.l_0.^2 .* Network0.E .* (r./Network0.l - 1) .* d ./ r;
f = sqrt(sum(f_s.^2,2))  .*  sign( r./Network0.l - 1);
PE = sum( Network0.l_0.^3 .* Network0.E .* (r./Network0.l - 1).^2 );

hfig = figure(149);  clf;  hfig.Position = [ 17         58        1844         930];
% tiles = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
% nexttile(1);
swimmer_axes(1) = axes('position',[ 0.041618       0.5871      0.62745      0.36228 ]);

swimmer_axes(2) = axes('position',[ 0.48609      0.57832      0.74588      0.40658]);


% 0.4897    0.4837    0.7110    0.4055

colors = repmat([0 0 0 ],Network0.n_nodes,1);  sizes = repmat(25,Network0.n_nodes,1);

% colors(Repulsion.is_inside == 1,:) = repmat([0 0 1 ],sum(Repulsion.is_inside == 1),1);
% sizes(Repulsion.is_inside == 1) = repmat(25,sum(Repulsion.is_inside == 1),1);

clear network_nodes network_links Mesh_faces Mesh_edges face_vertices0 edge_vertices0
for i = 1:length(swimmer_axes)
    axes(swimmer_axes(i));
    network_nodes(i) = scatter3(Network0.nodes(:,1),Network0.nodes(:,2),Network0.nodes(:,3),sizes,colors,'filled','o','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
    % pause
    hold on
    
    if plot_links
        network_links{i} = plot3([Network0.nodes(Network0.links(:,1),1) Network0.nodes(Network0.links(:,2),1)]',...
            [Network0.nodes(Network0.links(:,1),2) Network0.nodes(Network0.links(:,2),2)]',...
            [Network0.nodes(Network0.links(:,1),3) Network0.nodes(Network0.links(:,2),3)]',...
            '-','LineWidth',1.5,'Color',repmat(0.8,1,3));
    end
    
    
    [Mesh_faces{i},Mesh_edges{i}] = plot_mesh(Mesh0,2); Mesh_faces{i}(1).FaceAlpha = 1;
    Mesh_edges{i}(1).EdgeAlpha = 0.05;    Mesh_edges{i}(2).EdgeAlpha = 0.05;
    grid off
    
    
    xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
    
    
    hold off
    light
    
    axis equal
    %     view(-270,0);
    % view([0 90]);
    % view([-45.6 17]);
    % view([  -5.5940   19.6228]);
    
    xlim([-5.3 22]); ylim([-4.25 4.25]); zlim([-4.25 4.25]);
    
    
    % ax1_2 = copyobj(ax1,hfig);
    % ax1_2.View = [90 0];
    % ax1_2.Position = [0.4366    0.4837    0.7867    0.4754];
    
    
    if i == 1
%         title_handle = title(swimmer_axes(1),{"t = " + num2str(t) , "# offenders = " + num2str(sum(Repulsion.is_inside == 1))}); %,"PE = " + num2str(PE,2)
   title_handle = title(swimmer_axes(1),{"t = " + num2str(t) }); %,"PE = " + num2str(PE,2)
  
    end
    
    
end

swimmer_axes(1).View = [0 90];
swimmer_axes(2).View = [90 0]; % from in front
% swimmer_axes(2).View = [-90 0]; % from behind

[face_vertices0] = {Mesh_faces{1}.Vertices};  [edge_vertices0] = {Mesh_edges{1}.Vertices};
% store initial values, which are the same for different views of the same data, so don't need to store separate copies

%%
% figure(5); clf;

% https://www.mathworks.com/matlabcentral/answers/459385-multiple-y-axes-on-single-x-axis

% nexttile(2)
ax2 = axes('position',  [ 0.068765      0.37856     0.8875       0.1551],'FontSize',11);
yyaxis left
% velocity_plot = plot(stored_output.time(1:1),sqrt(sum(stored_output.derivatives(1:1,1:3).^2,2)),'-','LineWidth',2);
velocity_plot = plot(stored_output.time,sqrt(sum(stored_output.derivatives(:,1:3).^2,2)),'b-','LineWidth',1.5);
ylabel('|U| (\mum s^{-1})','color','b');
% xlabel('time (s)');
xlim([t tfinal]);
ax2.XTickLabel = [];
ax2.YColor = 'b';

ax2.XTickMode = 'manual'; ax2.YTickMode = 'manual';
% ax2.YLim = [min(ax2.YTick), max(ax2.YTick)];
% ax2.YLim = [11.5 14.5];
ax2.XLimMode = 'manual';
YTick = ax2.YTick;

hold on
ylims = ylim;
time_line(1) = plot([0 0], ylims,'--','color',repmat(0.5,1,3),'linewidth',1);


logical_inds = ~isnan(stored_output.link_breakage_time);
clear breakages
breakages{1} = plot( repmat(stored_output.link_breakage_time(logical_inds),1,2)' , repmat(ylims,sum(logical_inds),1)',':','linewidth',0.5,'Color',[0 0 0 0.25]);
xlim([t tfinal]);

yyaxis right
efficiency_plot = plot(stored_output.time,  ...
    (sum(stored_output.derivatives(:,1:3).^2 , 2)) ./ (stored_output.derivatives(:,7) .* ...
    assembly_input.Tail.motor_torque(stored_output.time) )   ,'r-','linewidth',1.5);
ylabel('|U|^2 / power','color','r');
xlim([t tfinal]);
ax2.YColor = 'r';
xlim([t tfinal]);
% ylim([7.3 14.2]);
% xlim([t tfinal]);


hold off
% yyaxis right
% PE_plot = plot(stored_output.time(1:1),PE(1),'-','linewidth',2);



ax2_1 = axes('Position',ax2.Position);

freq_plot = plot(stored_output.time,stored_output.derivatives(:,7),'-','linewidth',1.5);
freq_plot.Color = [0 0.65 0];
xlim([t tfinal]);
% ylim([-0.005 1.5E3]);
ax2_1.Color = 'none';


ax2_1.XTick = [];
% ax2.YLimMode = 'manual';
% ax2_1.YTick = linspace(ax2_1.YLim(1),ax2_1.YLim(2),length(YTick));
ax2_1.YTickLabel = strcat(ax2_1.YTickLabel, {'               '});


ylabel('tail frequency (rad/s)','color',[0 0.65 0]);
ax2_1.FontSize = 11;
ax2_1.YColor = [0 0.65 0];


ax3 = axes('position',  [  0.068765       0.20844           0.8875       0.1551 ],'FontSize',11);

PE_plot = plot(stored_output.time,stored_output.PE,'-','linewidth',2);
% ylim([-0.005 1.5E3]);
ylims = ylim;
xlim([t tfinal]);
ax3.XTickLabel = [];
hold on
time_line(2) = plot([0 0], ylims,'--','color',repmat(0.5,1,3),'linewidth',1);
breakages{2} = plot( repmat(stored_output.link_breakage_time(logical_inds),1,2)' , repmat(ylims,sum(logical_inds),1)',':','linewidth',0.5,'Color',[0 0 0 0.25]);

hold off
ylabel('PE');


ax4 = axes('position',  [ 0.068765        0.048219       0.8875       0.1551 ],'FontSize',11);
yyaxis left
% repulsion_force_plot = plot(stored_output.time,stored_output.repulsion.total_force_mag,'-','linewidth',2);
xlabel('time (s)');  xlim([0 tfinal]); ylims = ylim; ylim([-2 ylims(2)]);
ylabel('repulsive force');
hold on
ylims = ylim;
time_line(3) = plot([0 0], ylims,'--','color',repmat(0.5,1,3),'linewidth',1);
breakages{3} = plot( repmat(stored_output.link_breakage_time(logical_inds),1,2)' , repmat(ylims,sum(logical_inds),1)',':','linewidth',0.5,'Color',[0 0 0 0.25]);
hold off
yyaxis right
% repulsion_force_plot = plot(stored_output.time,stored_output.repulsion.total_torque_mag,'-','linewidth',2);
xlim([0 tfinal]);
ylim([-2 Inf]);
ylabel('repulsive torque');

drawnow


%% helps choose colormap limits for links


% N = 10;
% 
% f_all = NaN(Network.n_links,N);
% c = 0;
% for t = linspace(0,tfinal,N)  %  sol(end).x(end)
%     c = c + 1;
%     y = deval(sol, t);
%     
%     Network = Network0;  Network.E(stored_output.link_breakage_time <= t) = 0;
%     [Mesh, Network, mesh_node_parameters, y_swimmer] = update_Mesh_Network(Mesh0, Network,t, y, mesh_node_parameters,index_mapping, assembly_input);
%     
%     d = Network.nodes(Network.links(:,2),:)  - Network.nodes(Network.links(:,1),:) ;
%     r = sqrt(sum( d.^2 , 2 ) );
%     f_s = Network.l_0.^2 .* Network.E .* (r./Network.l - 1) .* d ./ r;
%     f = sqrt(sum(f_s.^2,2))  .*  sign( r./Network.l - 1);
%     f_all(:,c) = f;
% end
% figure(934)
% histogram(f_all(:))

%%
tic


save_vid = false;
if save_vid
    try, close(vidh); end
    %     vidh = VideoWriter('C:\Hull\swimmer network videos\swimmer base case','MPEG-4');
    vidh = VideoWriter('C:\Users\rudi\Desktop\RD\Tulane vids\link breakage body 0 tail 0 E 100 eta 5000 long network wider tube choppy maxstep 1E-3','MPEG-4');  % Motion JPEG AVI
    vidh.FrameRate = 30; %50; 30
    vidh.Quality = 95; %1-100
    resolution = '-r0';
    %     resolution = '-r180';  % doesn't seem to look much different, but
    %     doubles file size
    
    open(vidh);
end

clear links

  T = linspace(0,tfinal, 5); % 600   % 
% T = sol(1).x( sol(1).x <= Inf);
%  T = 0: 1*  1E-3 : tfinal;

time = T(  T >= 0 & T <= Inf);
% PE = NaN(1, length(time));


c= 0;
for t = time  %  sol(end).x(end)
    c = c+1;
    if t < 0.15
        ind = 1;
    else
        ind = 1;
    end
    
    y = deval(sol(ind), t);
    
    Network = Network0;  Network.E(stored_output.link_breakage_time <= t) = 0;
    [Mesh, Network, mesh_node_parameters, y_swimmer] = update_Mesh_Network(Mesh0, Network,t, y, mesh_node_parameters,index_mapping, assembly_input);
    
    
    %%
           active = true(Network.n_nodes,1); % active if at least one working link, passive if no active links (all E == 0)
        for i = 1:Network.n_nodes
            active(i) = any( Network.E(  Network.link_members{i}(1,:) ) ~= 0 ); % node is active if it has at least one link with nonzero stiffness
        end
        
%         active_network_nodes = Network.nodes(active,:);
        % any passive nodes default (and remain) as NaN / "not inside" (even though they actually may be inside)
        distance2mesh = NaN(Network.n_nodes,1);
        is_inside = false(Network.n_nodes,1);
        x = NaN(Network.n_nodes,3);
        
        [distance2mesh(active), x(active,:), xi_eta, mesh_index, element_index, is_inside(active)] = dist2mesh(Network.nodes(active,:), Mesh, index_mapping, assembly_input.repulsion.mindist2);
        % toc
        
        
        inside_threshold = 0.006; % body radius / 100
        % if inside further than this, stop simulation and put node back on surface
        
        distance2mesh(is_inside)
  
        
        if any(distance2mesh(is_inside) > 0.006 )
           bad_times(end+1) = t;
        end
        %%
    
    %     figure(351);
    %     histogram(f);
    
    
    %     dists = sqrt( sum((Network.nodes - Mesh(1).refpoints(:,1)').^2 , 2) ) - Metadata(1).geom.sphererad;
    %     insides = find(dists < 0);
    %     rel_dists = dists / Metadata(1).geom.sphererad;
    
    %     Repulsion = calc_repulsive_forces(Mesh, Network, assembly_input.repulsion);
    
    
    
    %     hfig = figure(142);
    %     axes(ax1);
    
    y_swimmer(7) = mod(y_swimmer(7) + 2*pi, 2*pi);  % shift periodic tail phase angle back to 0 - 2*pi
    
    refpoint = Mesh0(2).refpoints(:,1);
    % tail rotation
    rotmat = rotate_arbitrary_vector( Mesh0(assembly_input.Tail.submesh_index).orientation(:,1)  , y_swimmer(7));
    
    face_vertices = face_vertices0;  edge_vertices = edge_vertices0;
    
    face_vertices{2} = ( ( rotmat * (face_vertices0{2}' - refpoint) ) + refpoint )';
    edge_vertices{2} = ( ( rotmat * (edge_vertices0{2}' - refpoint) ) + refpoint )';
    % (body + tail) rotation
    rotmat = A_1_matrix(y_swimmer(4:6));
    for i = 1:2
        face_vertices{i} = ( ( rotmat * (face_vertices{i}' - refpoint) ) + refpoint )';
        edge_vertices{i} = ( ( rotmat * (edge_vertices{i}' - refpoint) ) + refpoint )';
    end
    
    for i = 1:2
        face_vertices{i} = face_vertices{i} + y_swimmer(1:3)';
        edge_vertices{i} = edge_vertices{i} + y_swimmer(1:3)';
    end
    
    nodes = reshape( y(Network0.n_links + 1 : Network0.n_links + Network0.n_nodes*3) , 3,[])';
    
    if plot_links
        links.x = num2cell([nodes(Network0.links(:,1),1) nodes(Network0.links(:,2),1)], 2);
        links.y =   num2cell([nodes(Network0.links(:,1),2) nodes(Network0.links(:,2),2)], 2);
        links.z =  num2cell([nodes(Network0.links(:,1),3) nodes(Network0.links(:,2),3)], 2);
        
        
        for ii = 1:2
            [network_links{ii}.XData] = deal(links.x{:}); [network_links{ii}.YData] = deal(links.y{:}); [network_links{ii}.ZData] = deal(links.z{:});
        end
        
        if color_links
            d = Network.nodes(Network.links(:,2),:)  - Network.nodes(Network.links(:,1),:) ;
            r = sqrt(sum( d.^2 , 2 ) );
            f_s = Network.l_0.^2 .* Network.E .* (r./Network.l - 1) .* d ./ r;
            f = sqrt(sum(f_s.^2,2))  .*  sign( r./Network.l - 1);
            %     PE(c) = sum( Network.l_0.^3 .* Network.E .* (r./Network.l - 1).^2 );
            
            [colors] = colordata(100,'whitejet',[-1 1],f); % was [-0.05 0.05]
            colors(:,4) = 0.3;
            
            links.color = num2cell(colors,2);
            
             links.visible = repmat( {"on"} , Network.n_links , 1);
        links.visible(t > stored_output.link_breakage_time) = {"off"};
        
            for ii = 1:2
                [network_links{ii}.Color] = deal( links.color{:});
                 [network_links{ii}.Visible] = deal( links.visible{:});
            end
            
        end
        
       
           
        
    end
    
    for ii = 1:2 % two views in two axes
        [Mesh_faces{ii}.Vertices] = face_vertices{:};  [Mesh_edges{ii}.Vertices] = edge_vertices{:};
        network_nodes(ii).XData = nodes(:,1)'; network_nodes(ii).YData = nodes(:,2)'; network_nodes(ii).ZData = nodes(:,3)';
    end
    %     tic
    %     links.x = [nodes(Network0.links(:,1),1) nodes(Network0.links(:,2),1)]';
    %     links.y =   [nodes(Network0.links(:,1),2) nodes(Network0.links(:,2),2)]';
    %     links.z =  [nodes(Network0.links(:,1),3) nodes(Network0.links(:,2),3)]';
    %
    %     for i = 1:Network0.n_links
    %         network_links(i).XData = links.x(:,i)'; network_links(i).YData = links.y(:,i)'; network_links(i).ZData = links.z(:,i)';
    %     end
    %     toc
    
    %     tic
    
    %     toc
    
    % filter = stored_output.time <= t;
    % velocity_plot.XData = stored_output.time(filter)';
    % velocity_plot.YData = sqrt(sum(stored_output.derivatives(filter,1:3).^2,2))';
    %
    % PE_plot.XData = time;
    % PE_plot.YData = PE;
    for iii = 1:3
        time_line(iii).XData = repmat(t,1,2);
    end
    
    [~,ind] = near(t,stored_output.time);
    any_inside = stored_output.trespassing.any_inside(ind);
    if any_inside
        str = "Network node(s) inside";
    else
        str = "";
    end
    
    title_handle.String = ({"t = " + num2str(t) , str});  % ,"PE = " + num2str(PE,2)
    if ~save_vid
        drawnow
    end
    %         pause
    
    if save_vid
        orig_mode = get(hfig, 'PaperPositionMode');
        
        set(hfig, 'PaperPositionMode', 'auto');
        cdata = print(hfig,'-RGBImage','-opengl',resolution);
        % Restore figure to original state
        set(hfig, 'PaperPositionMode', orig_mode); % end
        
        %For the "OpenGL" renderer you can write a similar code. This technique will not work for the "painters" renderer.
        
        %Next, replace the use of GETFRAME from your code with IM2FRAME as follows:
        
        F = im2frame(cdata);
        
        writeVideo(vidh,F);
    end
    
    
end

if save_vid
    close(vidh);
end

toc