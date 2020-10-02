load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5_AR2_0.5_tail_radius_0.03101752454497_amp_0.62734375_lambda_4.390625_nlambda_1.35703125_motorBC_torque_dump.mat');
load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5_AR2_0.5_tail_radius_0.03101752454497_amp_0.62734375_lambda_4.390625_nlambda_1.35703125_motorBC_torque_timestepping.mat');

% load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5_AR2_0.5_tail_radius_0.03101752454497_amp_0.62734375_lambda_4.390625_nlambda_1.3565625_motorBC_torque_dump.mat');
% load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5_AR2_0.5_tail_radius_0.03101752454497_amp_0.62734375_lambda_4.390625_nlambda_1.3565625_motorBC_torque_timestepping.mat');

% load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5_AR2_0.5_tail_radius_0.03101752454497_amp_0.627_lambda_4.39_nlambda_1.357_motorBC_torque_dump.mat');
% load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5_AR2_0.5_tail_radius_0.03101752454497_amp_0.627_lambda_4.39_nlambda_1.357_motorBC_torque_timestepping.mat');

% load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5_AR2_0.7_tail_radius_0.03101752454497_amp_0.627_lambda_4.39_nlambda_1.357_motorBC_torque_dump.mat');
% load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5_AR2_0.7_tail_radius_0.03101752454497_amp_0.627_lambda_4.39_nlambda_1.357_motorBC_torque_timestepping.mat');

% load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5.75_AR2_0.65_tail_radius_0.03101752454497_amp_0.634_lambda_4.32_nlambda_1.345_motorBC_torque_dump.mat');
% load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5.75_AR2_0.65_tail_radius_0.03101752454497_amp_0.634_lambda_4.32_nlambda_1.345_motorBC_torque_timestepping.mat');

% load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5.75_AR2_0.65_tail_radius_0.03101752454497_amp_0.6335766697050216_lambda_4.321434829388381_nlambda_1.344662931002491_motorBC_torque_dump.mat');
% load('C:\Hull\methods_fig_dumps\curved_rod_AR1_5.75_AR2_0.65_tail_radius_0.03101752454497_amp_0.6335766697050216_lambda_4.321434829388381_nlambda_1.344662931002491_motorBC_torque_timestepping.mat');
fig_type = 'path';
% fig_type = 'diffusivities';

tracked_pt = timestepping_solution.refpoint; % center of body
% tracked_pt = timestepping_solution.y(:,1:3);   % more or less ass end

 [shift,angle,rotvec] = align_path(fits,tracked_pt(1,:));
rotmat_align = rotate_arbitrary_vector(rotvec,angle);
aligned_path = tracked_pt' + repmat(shift',1,size(tracked_pt,1));
aligned_path = (rotmat_align * aligned_path)';
switch fig_type
    case 'path'
do_inset = true;
do_path = true;
    case 'diffusivities'
        do_inset = false;
        do_path = false;
end
ind = 1015;  % timestepping ind to draw body at and path up till
start_ind = 140;


figure(456);
clf;
set(gcf,'Position',[    350         -71        1439        1041]);
[Mesh2 , rotmat_time] = move_Mesh(Mesh, timestepping_solution.y(ind,1:6));
Mesh2 = shiftMesh(Mesh2, shift);
[Mesh2 , rotmat_align] = rotateMesh(Mesh2,angle,rotvec);
[s,e] = plot_mesh(Mesh2,[3 4]);
e(2).EdgeAlpha = 0;
switch fig_type
    case 'path'
        e(1).EdgeAlpha = 0.3;
    case 'diffusivities'
        e(1).EdgeAlpha = 0.3;
        s(1).FaceAlpha = 0.7;
         s(2).FaceAlpha = 0.25;
end

% set(gca,'View',[-30 24.4]);
hold on
if do_path 
% hp = plot3(timestepping_solution.refpoint(1:ind,1),timestepping_solution.refpoint(1:ind,2),timestepping_solution.refpoint(1:ind,3),'k-','linewidth',2);
hpp = plot3(aligned_path(start_ind:ind,1),aligned_path(start_ind:ind,2),aligned_path(start_ind:ind,3),'k-','linewidth',1.5);
end

transform = @(vec) rotmat_align * (  (rotmat_time * vec) + timestepping_solution.y(ind,1:3)' + shift'  ); %contains all necessary rotations, translations accounting for path alignment and moving to time represented by ind
rotmat_tot = rotmat_align * rotmat_time;

switch fig_type
    case 'diffusivities'
%% center of diffusion

center0 = D.center + Mesh(2).refpoints;  % for whatever reason, D.center defined relative to tail refpoint when tail is present, body refpoint when tail isn't present....
Cd_transformed = transform(center0);
% cd = plot3(Cd_transformed(1),Cd_transformed(2),Cd_transformed(3),'ko','markerfacecolor','k','markersize',10);
try, delete(cd), end;
cd = drawSphere([Cd_transformed' 0.025],'nPhi',50,'nTheta',50);
cd.EdgeColor = 'none';  cd.FaceLighting = 'none';  cd.FaceColor = 'k';

main_ax = gca;

D2.rotation.diffusivity = diag(repmat(0.167,1,3));


%% principle axes of rotational diffusion
try delete(q); end
clear q

sf = [25 14 14];  maxheadsizes = [1*0.035 0.8*0.035 0.8*0.035];   colors_lines = [repmat(0.,1,3);  repmat(0.,2,3)];
styles = {'--','-','-'};
widths = 0.7*[0.02 0.015 0.015];
axes_transformed = rotmat_tot * D.rotation.axes;
fudge = [1.52 1 1];
%   arrow_ends = repmat(Cd_transformed,1,3) + repmat( sf*diag(D.rotation.diffusivity),1,3) .* axes_transformed;
      for qi = 1:3
q(qi) = plot3( [Cd_transformed(1) - fudge(qi)*sf(qi)*D2.rotation.diffusivity(qi,qi)*axes_transformed(1,qi) ; Cd_transformed(1) + sf(qi)*D2.rotation.diffusivity(qi,qi)*axes_transformed(1,qi) ] , ...
    [Cd_transformed(2) - fudge(qi)*sf(qi)*D2.rotation.diffusivity(qi,qi)*axes_transformed(2,qi) ; Cd_transformed(2) + sf(qi)*D2.rotation.diffusivity(qi,qi)*axes_transformed(2,qi) ] , ...
    [Cd_transformed(3) - fudge(qi)*sf(qi)*D2.rotation.diffusivity(qi,qi)*axes_transformed(3,qi) ; Cd_transformed(3) + sf(qi)*D2.rotation.diffusivity(qi,qi)*axes_transformed(3,qi) ] );
set(q(qi),'Color',[colors_lines(qi,:) 0.15],'LineStyle',styles{qi},'LineWidth',1);
      end


% for qi = 1:3
% %     q(qi) = quiver3(Cd_transformed(1),Cd_transformed(2),Cd_transformed(3),...
% %         sf*D.rotation.diffusivity(qi,qi)*axes_transformed(1,qi),...
% %         sf*D.rotation.diffusivity(qi,qi)*axes_transformed(2,qi),...
% %         sf*D.rotation.diffusivity(qi,qi)*axes_transformed(3,qi),...
% %         0,'color',colors(qi),'linewidth',3,'maxheadsize',maxheadsizes(qi));
%     
%         q(qi) = mArrow3( [Cd_transformed(1),Cd_transformed(2),Cd_transformed(3)],...
%         [Cd_transformed(1) + sf*D2.rotation.diffusivity(qi,qi)*axes_transformed(1,qi),...
%         Cd_transformed(2) + sf*D2.rotation.diffusivity(qi,qi)*axes_transformed(2,qi),...
%         Cd_transformed(3) + sf*D2.rotation.diffusivity(qi,qi)*axes_transformed(3,qi)],...
%         'color',colors(qi,:),'stemWidth',widths(qi),'tipWidth',maxheadsizes(qi));
% end
% 
% for qii = 1:3
% %     q(qi) = quiver3(Cd_transformed(1),Cd_transformed(2),Cd_transformed(3),...
% %         sf*D.rotation.diffusivity(qi,qi)*axes_transformed(1,qi),...
% %         sf*D.rotation.diffusivity(qi,qi)*axes_transformed(2,qi),...
% %         sf*D.rotation.diffusivity(qi,qi)*axes_transformed(3,qi),...
% %         0,'color',colors(qi),'linewidth',3,'maxheadsize',maxheadsizes(qi));
%     
%         q(qii) = mArrow3( [Cd_transformed(1),Cd_transformed(2),Cd_transformed(3)],...
%         [Cd_transformed(1) - sf*D2.rotation.diffusivity(qii,qii)*axes_transformed(1,qii),...
%         Cd_transformed(2) - sf*D2.rotation.diffusivity(qii,qii)*axes_transformed(2,qii),...
%         Cd_transformed(3) - sf*D2.rotation.diffusivity(qii,qii)*axes_transformed(3,qii)],...
%         'color',colors(qii,:),'stemWidth',widths(qii),'tipWidth',maxheadsizes(qii));
% end
% set(gca,'view',[-1.898437500000000   4.729492187500000]);
% set(gca,'CameraPosition',[ -41.645793144791604 -37.769559708763801   3.084233709554325]);
% set(gca,'CameraTarget',[  8.339591767865066   0.171681052500000   0.015260851427763]);
% set(gca,'CameraViewAngle',   5.450558963549248);
axes(main_ax);


%%
axes(main_ax);
try, delete(ch); end
try, delete(cha); end
try, delete(Cones); end; try, delete(EndPlate1); end; try, delete(EndPlate2); end;
try, delete(Cones_2); end; try, delete(EndPlate1_2); end; try, delete(EndPlate2_2); end;


hold on
theta= {2*pi*0.005 :0.02:   2*pi * 0.995;  2*pi*0.02 :0.02:   2*pi * 0.98;  2*pi*0.02 :0.02:   2*pi * 0.98;};
% theta = 0:0.02:2*pi;
% radii = [0.6 0.227 0.234];
for qi = 1:3
v=null(axes_transformed(:,qi)');
radius = 12 * D.rotation.diffusivity(qi,qi);
% radius = radii(qi);
%  radius = 0.165;  %radius = 0.125;
 frac = [0.9 0.85 0.7];  frac = [0 0 0];
 flips = [1 1 1];  rot_angles = [60*pi/180 105*pi/180 228*pi/180];
 cone_widths = [0.025 0.02 0.02];
 cone_inds = [4 1; 20 1; 20 1];
 colors = [0 0 0; 0 0 0; 0 0 0];

 center = [Cd_transformed(1) + flips(qi)*frac(qi)*sf(qi)*D.rotation.diffusivity(qi,qi)*axes_transformed(1,qi) ; ...
           Cd_transformed(2) + flips(qi)*frac(qi)*sf(qi)*D.rotation.diffusivity(qi,qi)*axes_transformed(2,qi) ; ...
           Cd_transformed(3) + flips(qi)*frac(qi)*sf(qi)*D.rotation.diffusivity(qi,qi)*axes_transformed(3,qi) ];
       
%         rotmat = rotate_arbitrary_vector( axes_transformed(:,qi), rot_angles(qi));
points= repmat(center,1,size(theta{qi},2))+radius*(v(:,1)*cos(theta{qi})+v(:,2)*sin(theta{qi}));
points = ( rotate_arbitrary_axis(points', Cd_transformed, axes_transformed(:,qi), rot_angles(qi)) )';
ch(qi) = plot3(points(1,4:end-4),points(2,4:end-4),points(3,4:end-4),'linewidth',0.75,'LineStyle',styles{qi},'Color',colors(qi,:));

%  cha(qi) = mArrow3( points(:,1),points(:,2), 'color','k','stemWidth',10,'tipWidth',20);
[Cones(qi),EndPlate1(qi),EndPlate2(qi)] = Cone(points(:,cone_inds(qi,1)),points(:,cone_inds(qi,2)),[cone_widths(qi) 0],300,colors(qi,:),1,true,false);
[Cones_2(qi),EndPlate1_2(qi),EndPlate2_2(qi)] = Cone(points(:,end-cone_inds(qi,1)+1),points(:,end-cone_inds(qi,2)+1),[cone_widths(qi) 0],300,colors(qi,:),1,true,false);
end



%   q(qi) = mArrow3( [Cd_transformed(1),Cd_transformed(2),Cd_transformed(3)],...
%         [Cd_transformed(1) + sf*D.rotation.diffusivity(qi,qi)*axes_transformed(1,qi),...
%         Cd_transformed(2) + sf*D.rotation.diffusivity(qi,qi)*axes_transformed(2,qi),...
%         Cd_transformed(3) + sf*D.rotation.diffusivity(qi,qi)*axes_transformed(3,qi)],...
%         'color',colors(qi,:),'stemWidth',widths(qi),'tipWidth',maxheadsizes(qi));
    
    
% ch(1) = drawCircleArc3d([Cd_transformed'  3     45 0        0   0 350]);
% ch(2) = drawCircleArc3d([Cd_transformed'  3/9     90 0        0   0 350]);
% ch(3) = drawCircleArc3d([Cd_transformed'  3/9     0 90        0   0 350]);
%%
% function varargout = drawCircleArc3d(arc, varargin)
%DRAWCIRCLEARC3D Draw a 3D circle arc
%
%   drawCircleArc3d([XC YC ZC R THETA PHI PSI START EXTENT])
%   [XC YC ZC]  : coordinate of arc center
%   R           : arc radius
%   [THETA PHI] : orientation of arc normal, in degrees (theta: 0->180).
%   PSI         : roll of arc (rotation of circle origin)
%   START       : starting angle of arc, from arc origin, in degrees
%   EXTENT      : extent of circle arc, in degrees (can be negative)
%   




end



set(gca,'view',[-31.875390625 15.928515625]);
set(gca,'CameraPosition',[ -23.828303617479669 -47.212434274339280  16.019837180219831]);
set(gca,'CameraTarget',[   5.637032765231945   0.171681052500000   0.095013656500000]);
set(gca,'CameraViewAngle',   6.25954081584161);



% set(axes1,'CameraViewAngle',6.25954081584161,'DataAspectRatio',[1 1 1],...
%     'PlotBoxAspectRatio',[5.82576212474715 1.01165408776581 1]);


axis equal
axis tight
 axis off
 l = light('position',[-41.645793144791604 -37.769559708763801   3.084233709554325]);

%  uistack(hpp,'bottom');  doesn't do anything, I guess 3D position is
%  different than stacking

%% inset of path viewed end-on
if do_inset
   axes3 = axes('Parent',gcf,...
        'Position',[(0.080491095686286 -0.08)   (0.477700373831766 -0.02)  0.235468904313713   0.105653644859812]); 
    axes3.Units = 'centimeters';  axes3.Position = [8.6503 + 1 15.123 2.53 2.53];
    
    hp_inset = plot3(aligned_path(1:ind,1),aligned_path(1:ind,2),aligned_path(1:ind,3),'k-','linewidth',1.5);
    axis equal;  axis tight;  axis off;
    view([-90 0]);
    
   ann = annotation(gcf,'rectangle',...
        [0.2272 0.56761 0.095 0.095] ,...
        'LineWidth',3);
    ann.Units = 'centimeters';  ann.Position = [8.4122 + 1       14.938         2.95         2.95];
end

view(main_ax,[-65 5.555]);

return
% this is the wrong function, should be export_fig.....
% exportfig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\curved rod path and principle axes.png','format','png','color','rgb','resolution',800)
%%
export_fig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\curved rod swimming path.png','-r800','-transparent')
%%
export_fig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\curved rod diffusivities.png','-r800','-transparent')