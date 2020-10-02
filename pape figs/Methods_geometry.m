dumpfolder = 'C:\Hull\select_dumps\';
File = 'curved_rod_AR1_5_AR2_0.5_tail_radius_0.03101752454497_amp_0.6203184964759591_lambda_4.239738315401703_nlambda_1.337042130119591_motorBC_torque_dump.mat';

load([dumpfolder,File],'Metadata','Mesh');
% Mesh = move_Mesh(Mesh,[0 0 0 pi/2 * 0.65 0 0 0]');
% rotangle = [pi/2 * 0.65 0 0];

body_only = true;

do_text = false;  % do text labels?

Mesh = move_Mesh(Mesh,[0 0 0    0 0 pi/2 * 0.65]');
rotangle = [0 0 pi/2 * 0.65];

figure(926)
set(gcf,'Position',[ 109         -39        1715        1021]);
% set(gcf,'Renderer','painters');
clf
set(gcf,'Color','w');

if body_only
    [s,e] = plot_mesh(Mesh(1));
else
    [s,e] = plot_mesh(Mesh);
end
hold on
set(gca,'view',[   -47    12]);

% set(gcf,'position',[          338         -32        1452        1018]);
%   set(gcf,'Position',[  2053         -60        1452         999]);
axis tight
axis equal
light
lighting phong

set(s,'ambientStrength',0.7);
set(e,'edgealpha',0.075);
set(s,'facealpha',0.4);
set(s,'facecolor',repmat(0.7,1,3));



r = Metadata(1).geom.radius_curv;
centx = 0;  centy = r;
%theta = 0 : (2 * pi / 1000) : (2 * pi);
theta = pi*0.4 : (2 * pi / 1000) : (1.38 * pi);
pline_x = r * cos(theta) + centx;
pline_y = r * sin(theta) + centy;
pline_z = zeros(size(pline_x)) + 0;
circpts = [pline_x; pline_y; pline_z];
rotmat = rotation_matrix('z',rotangle(1)) * rotation_matrix('y',rotangle(2)) * rotation_matrix('x',rotangle(3));
circpts1 = rotmat * circpts;
% circpts1 =  circpts;

circ1 = plot3(circpts1(1,:),circpts1(2,:),circpts1(3,:),':k','linewidth',2);

theta = (1.38 * pi) : (2 * pi / 1000) : 2*pi+pi*0.4;
pline_x = r * cos(theta) + centx;
pline_y = r * sin(theta) + centy;
pline_z = zeros(size(pline_x)) + 0;
circpts = [pline_x; pline_y; pline_z];
rotmat = rotation_matrix('z',rotangle(1)) * rotation_matrix('y',rotangle(2)) * rotation_matrix('x',rotangle(3));
circpts2 = rotmat * circpts;
%   circpts2 = circpts;

circ2 = plot3(circpts2(1,:),circpts2(2,:),circpts2(3,:),'-k','linewidth',2);

arrow1 = rotmat * [ [0 r 0]' [r r 0]' ];
% arrow1 = [ [0 r 0]' [r r 0]' ];
% ah1 = arrow(arrow1(:,1),arrow1(:,2),'ends','both','length',25,'width',2);
% ah1a = arrow3(arrow1(:,1)',arrow1(:,2)');
stemWidth = 0.005;  tipWidth = 0.04;
ah1a = mArrow3(arrow1(:,1)',arrow1(:,2)','stemWidth',stemWidth,'tipWidth',tipWidth);
ah1b = mArrow3(arrow1(:,2)',arrow1(:,1)','stemWidth',stemWidth,'tipWidth',tipWidth);
%ah1b = arrow3(arrow1(:,2)',arrow1(:,1)','linewidth',1.5,'width',1E-2);
%arrow1 = rotmat * [ [0 2*r 0]' [0 r 0]' ];

theta = 0 : (2 * pi / 1000) : (2 * pi);
pline_x = r * cos(theta) + centx;
pline_y = r * sin(theta) + centy;
pline_z = zeros(size(pline_x)) + 0;
circpts = [pline_x; pline_y; pline_z];
rotmat = rotation_matrix('z',rotangle(1)) * rotation_matrix('y',rotangle(2)) * rotation_matrix('x',rotangle(3));
circpts = rotmat * circpts;

% ah2 = arrow(circpts(:,212),circpts(:,213),'length',25);
% ah3 = arrow(circpts(:,687),circpts(:,686),'length',25);
stemWidth = 0.005;  tipWidth = 0.04;
ah2 = mArrow3(circpts(:,180)',circpts(:,213)','stemWidth',stemWidth,'tipWidth',tipWidth);
ah3 = mArrow3(circpts(:,719)',circpts(:,686)','stemWidth',stemWidth,'tipWidth',tipWidth);


if ~body_only
    % xdirh = arrow(rotmat*Mesh(2).refpoints,rotmat*(Mesh(2).refpoints + [1 0 0]'),'length',25,'width',2,'type','line','linestyle','--');
    % ydirh = arrow(rotmat*Mesh(2).refpoints,rotmat*(Mesh(2).refpoints + [0 1 0]'),'length',25,'type','line','linestyle','--');
    % zdirh = arrow(rotmat*Mesh(2).refpoints,rotmat*(Mesh(2).refpoints + [0 0 1]'),'length',25,'width',2,'type','line','linestyle','--');
    stemWidth = 0.005;  tipWidth = 0.02;
    xdirh = mArrow3(rotmat*Mesh(2).refpoints,rotmat*(Mesh(2).refpoints + [1.5*1.3 0 0]'),'stemWidth',stemWidth,'tipWidth',tipWidth);
    ydirh = mArrow3(rotmat*Mesh(2).refpoints,rotmat*(Mesh(2).refpoints + [0 1*1.3 0]'),'stemWidth',stemWidth,'tipWidth',tipWidth);
    zdirh = mArrow3(rotmat*Mesh(2).refpoints,rotmat*(Mesh(2).refpoints + [0 0 1.5*1.3]'),'stemWidth',stemWidth,'tipWidth',tipWidth);
    color_frac = 0.8;
    set(xdirh,'facecolor',repmat(color_frac,1,3)); set(ydirh,'facecolor',repmat(color_frac,1,3)); set(zdirh,'facecolor',repmat(color_frac,1,3));
    
    temp = rotmat*Mesh(2).refpoints; clear line;
    tl = line([temp(1) -6.5],[temp(2) 0],[temp(3) 0],'linestyle','--','color',repmat(color_frac,1,3));
    
    stemWidth = 0.005;  tipWidth = 0.04;
    lam1 = mArrow3([-1.37 -0.25 0]',[-1.37 - Metadata(2).geom.lambda -0.25 0]','stemWidth',stemWidth,'tipWidth',tipWidth);
    lam2 = mArrow3([-1.37 - Metadata(2).geom.lambda -0.25 0]',[-1.37 -0.25 0]','stemWidth',stemWidth,'tipWidth',tipWidth);
    
    stemWidth = 0.005;  tipWidth = 0.03;
    amp1 = mArrow3([-1.37 0 0]',[-1.37 0.188 0.5031]','stemWidth',stemWidth,'tipWidth',tipWidth);
    amp2 = mArrow3([-1.37 0.188 0.5031]',[-1.37 0 0]','stemWidth',stemWidth,'tipWidth',tipWidth);
end

stemWidth = 0.005;  tipWidth = 0.03;
d1 = mArrow3([0.5217 0.9763 1.259]',[0.9295 0.92 1.758]','stemWidth',stemWidth,'tipWidth',tipWidth);
d2 = mArrow3([0.9295 0.92 1.758]',[0.5217 0.9763 1.259]','stemWidth',stemWidth,'tipWidth',tipWidth);

if ~body_only
    X0 = plot3(Mesh(1).refpoints(1),Mesh(1).refpoints(2),Mesh(1).refpoints(3),'o','markerfacecolor',repmat(0.3,1,3),'MarkerEdgeColor',repmat(0.3,1,3),'MarkerSize',10);
end

coords = [  0.55695      0.47863     0.020802     0.027484;...  a
    0.50178       0.3646     0.042605     0.033827;... lambda
    0.6571      0.62287     0.041151     0.029598;...
    0.69674      0.69711     0.039698     0.035941;...  l
    0.76275      0.49719     0.023709     0.021142;...  x
    0.75921      0.56927     0.045512     0.027484;...  y
    0.59358      0.61723     0.022256     0.027484;...  z
    0.66251      0.77456     0.039698      0.02537;...  d
    0.71781      0.61442     0.045512     0.027484;...  x0
    ];

coords_body = [   0.450871635829806   0.614748677118797   0.915399177675431;
    0.736557624677545   0.494724988060362   1.316331963546851;
    0.548 0.677 1.64];

if do_text
    
    texts = {'$$ a $$', '$$ \lambda $$','$$ r $$','$$ \ell $$','$$ x $$','$$ y $$','$$ z $$','$$ d $$', '{\boldmath$X_0$}'};
    texts_body = {'$$ r $$','$$ \ell $$','$$ d $$'};
    
    if body_only
        fontsize_text = 34;
    else
        fontsize_text = 24;
    end
    
    if ~body_only
        for c = 1:length(texts)
            % Create textbox
            annotation(gcf,'textbox',...
                coords(c,:),...
                'String',{texts{c}},...
                'FitBoxToText','off','linestyle','none','Interpreter','latex','fontsize',fontsize_text);
        end
    else
        for c = 1:length(texts_body)
            % Create textbox
            text(coords_body(c,1),coords_body(c,2),coords_body(c,3),texts_body{c},...
                'Interpreter','latex','fontsize',fontsize_text,'fontweight','bold');
        end
        
        %     SFdefs = text( -1.1443 ,      2.4353  ,     1.2612,{'SF_1 = \itl / \itd','\rmSF_2 = \itl / ( 2\pi\itr )'},'fontsize',18,'edgecolor','none','linewidth',1.5,'interpreter','tex');
        SFdefs = text(   -0.213047321043767  , 1.404318987290871 ,  1.157451610181521,{'$$ \mathcal{L} = \ell / d \smallskip $$','$$ \mathcal{K} = \ell / (2\pi r) $$'},'fontsize',fontsize_text,'edgecolor','none','linewidth',1.5,'interpreter','latex');
        % SFdefs = text( -1.1443 ,      2.4353  ,     1.2612,{'shat $$ \mathcal{L} $$','crat'},'fontsize',18,'edgecolor','none','linewidth',1.5,'interpreter','latex');
    end
    
end
%{'SF_1 = l / d','SF_2 = l / (2\pir)','$\e^{\frac{1}{2x}}$'}
axis off

if ~body_only
    %% inset body mesh
    
    axes2 = axes('Parent',gcf,...
        'Position',[0.381931407967041 0.454507998868006 0.109804129223041 0.293139497985231]);
    
    [s,e] = plot_mesh(Mesh(1));
    hold on
    set(gca,'view',[   -47    12]);
    axis tight
    axis equal
    light
    lighting phong
    
    set(s,'ambientStrength',0.7);
    set(e,'edgealpha',1);
    set(s,'facealpha',1);
    set(s,'facecolor',repmat(0.7,1,3));
    axis off
    set(gca,'CameraPosition',[-10.899 -10.618 3.717],'CameraTarget',[0.67762 0.1774 0.35241]);
    
    %% inset tail mesh
    
    axes3 = axes('Parent',gcf,...
        'Position',[0.11134 0.51742 0.25462 0.18042]);
    
    [s,e] = plot_mesh(Mesh(2));
    hold on
    set(gca,'view',[   -47    12]);
    axis tight
    axis equal
    xlim([-6.1 -5.85])
    light
    lighting phong
    
    set(s,'ambientStrength',0.7);
    set(e,'edgealpha',1);
    set(s,'facealpha',1);
    set(s,'facecolor',repmat(0.7,1,3));
    axis off
    set(gca,'CameraPosition',[-7.385 -0.70457 0.71916],'CameraTarget',[-6.0381 0.55142 0.32771]);
    
    % inset rectangle outline
    annotation(gcf,'rectangle',...
        [0.21006 0.53445 0.25964 0.19291],...
        'LineWidth',3);
end

% for some reason export_fig fucks up the inset tail mesh, so have been
% using GUI saveas...
% export_fig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\survey SI.png','-r600','-nocrop')

export_fig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\methods notext.png','-r350','-p0.01')