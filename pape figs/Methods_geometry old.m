dumpfolder = 'E:\Hull\select_dumps\';
File = 'curved_rod_AR1_5_AR2_0.5_tail_radius_0.03101752454497_amp_0.6203184964759591_lambda_4.239738315401703_nlambda_1.337042130119591_motorBC_torque_dump.mat';

load([dumpfolder,File],'Metadata','Mesh');
Mesh = move_Mesh(Mesh,[0 0 0 pi/2 * 0.65 0 0 0]');
rotangle = [pi/2 * 0.65 0 0];

figure(923)
clf
[s,e] = plot_mesh(Mesh);
hold on
 set(gca,'view',[   -47    12]);
 
 set(gcf,'position',[          338         -32        1452        1018]);
 axis tight
 axis equal
    light
    lighting phong
    
    set(s,'ambientStrength',0.8)
   
    
    r = Metadata(1).geom.radius_curv;
centx = 0;  centy = r;
%theta = 0 : (2 * pi / 1000) : (2 * pi);
theta = pi*0.4 : (2 * pi / 1000) : (1.38 * pi);
pline_x = r * cos(theta) + centx;
pline_y = r * sin(theta) + centy;
pline_z = zeros(size(pline_x)) + 0;
circpts = [pline_x; pline_y; pline_z];
  rotmat = rotation_matrix('z',rotangle(3)) * rotation_matrix('y',rotangle(2)) * rotation_matrix('x',rotangle(1));
  circpts1 = rotmat * circpts;
% circpts1 =  circpts;

circ1 = plot3(circpts1(1,:),circpts1(2,:),circpts1(3,:),':k','linewidth',2);

theta = (1.38 * pi) : (2 * pi / 1000) : 2*pi+pi*0.4;
pline_x = r * cos(theta) + centx;
pline_y = r * sin(theta) + centy;
pline_z = zeros(size(pline_x)) + 0;
circpts = [pline_x; pline_y; pline_z];
  rotmat = rotation_matrix('z',rotangle(3)) * rotation_matrix('y',rotangle(2)) * rotation_matrix('x',rotangle(1));
  circpts2 = rotmat * circpts;
%   circpts2 = circpts;

circ2 = plot3(circpts2(1,:),circpts2(2,:),circpts2(3,:),'-k','linewidth',2);

arrow1 = rotmat * [ [0 r 0]' [r r 0]' ];
% arrow1 = [ [0 r 0]' [r r 0]' ];
ah1 = arrow(arrow1(:,1),arrow1(:,2),'ends','both','length',25,'width',2);

%arrow1 = rotmat * [ [0 2*r 0]' [0 r 0]' ];

theta = 0 : (2 * pi / 1000) : (2 * pi);
pline_x = r * cos(theta) + centx;
pline_y = r * sin(theta) + centy;
pline_z = zeros(size(pline_x)) + 0;
circpts = [pline_x; pline_y; pline_z];
  rotmat = rotation_matrix('z',rotangle(3)) * rotation_matrix('y',rotangle(2)) * rotation_matrix('x',rotangle(1));
  circpts = rotmat * circpts;
ah2 = arrow(circpts(:,212),circpts(:,213),'length',25);

ah3 = arrow(circpts(:,687),circpts(:,686),'length',25);
