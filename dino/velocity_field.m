
field_input.accuracy.integration_tol.field_vel.abstol = 1E-7;
field_input.accuracy.integration_tol.field_vel.abstol = 1E-6;
%field_input.accuracy.integration_tol.field_vel.abstol = 1E-5;

field_input.accuracy.integration_tol.field_vel.reltol = 0;
field_input.accuracy.integration_tol.field_vel.maxevals = Inf;
field_input.accuracy.eps2 = input.accuracy.eps2;
field_input.constants.multfactor = input.constants.multfactor;
field_input.constants.mu = input.constants.mu;

%integrand_constants = struct('eps2',assembly_input.accuracy.eps2,'x_col',Mesh(i_mesh_vert).verts(local_vert,:)');


% load('C:\Users\rudi\Desktop\RD\dino_dumps_hairs\dino_parallel2_3_all_tailed3_body-transverse-tail-wingtip_dump.mat');

% load('C:\Users\rudi\Desktop\RD\dino_dumps_hairs\dino_parallel2_3_all_tailed3_wingtip_dump.mat');
% load('C:\Users\rudi\Desktop\RD\dino_dumps_hairs\dino_parallel2_3_all_tailed3_transverse_dump.mat');
% load('C:\Users\rudi\Desktop\RD\dino_dumps_hairs\dino_parallel2_3_all_tailed3_body-transverse-tail-wingtip_dump.mat');
load('C:\Users\rudi\Desktop\RD\dino_dumps_hairs\dino_refined_transverse_body-transverse-tail-wingtip_dump.mat');

% has Solutions struct
Solutions_fields = setdiff( fieldnames(Mesh_files) , {'phase','time','metadata'});  %order of submeshes as far as Solutions (e.g. traction f) is concerned
% hmm, it actually appears that the order of Solutions_fields is the same
% as the order of Solutions.rand_inds, so probably don't need to do
% anything to match the orders

% x = linspace(-2,3,10); % x = 1.6;
% y = linspace(-1.5,1.5,10);  
% z = linspace(-1.5,1.5,10); 
phase_ind = 1;  
phase = Solutions.phase(phase_ind);

 [Mesh, Metadata, matrix_props, Metadata_rand_inds] = load_dino_mesh(phase, Mesh_files, input, Solutions.rand_inds{phase_ind});

%%
Xlist = [];  Ylist = [];  Zlist = [];

pts_per_slice = 200;
pts_per_slice = 100;

Lx = 40;
Ly = 35;
Lz = 35;

% values at right are for current geom

% vertical planes through transverse

x = 7; %bottom of groove  % 6
y = linspace(-Ly,Ly,pts_per_slice);
z = linspace(-Lz,Lz,pts_per_slice);
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];  

x = 8.75; %middle of groove  % 9
y = linspace(-Ly,Ly,pts_per_slice);
z = linspace(-Lz,Lz,pts_per_slice);
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];  

x = 10.5; %top of groove % 12.5
y = linspace(-Ly,Ly,pts_per_slice);
z = linspace(-Lz,Lz,pts_per_slice);
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];  


% horizontal planes parallel to tail beating plane

z = 0;  %middle of body, below tail % 0
x = linspace(-Lx,Lx,pts_per_slice);
y = linspace(-Ly,Ly,pts_per_slice);
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];  

z = 15;  %bottom of tail % 10
x = linspace(-Lx,Lx,pts_per_slice);
y = linspace(-Ly,Ly,pts_per_slice);
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];  

z = 15.25;  %middle of tail  % 15
x = linspace(-Lx,Lx,pts_per_slice);
y = linspace(-Ly,Ly,pts_per_slice);
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];  

z = 15.5;  %top of tail % 20
x = linspace(-Lx,Lx,pts_per_slice);
y = linspace(-Ly,Ly,pts_per_slice);
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];  

% vertical, along tail axis, down through body

y = 0;  %through center of body, vertically, through tail % 0
x = linspace(-Lx,Lx,pts_per_slice);
z = linspace(-Lz,Lz,pts_per_slice);
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];  

y = -2;  %near center of body, vertically, slightly left of tail  % -2
x = linspace(-Lx,Lx,pts_per_slice);
z = linspace(-Lz,Lz,pts_per_slice);
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];  


y = 2;  %near center of body, vertically, slightly right of tail  % 2
x = linspace(-Lx,Lx,pts_per_slice);
z = linspace(-Lz,Lz,pts_per_slice);
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];  

%% whole volume field calc

Xlist = [];  Ylist = [];  Zlist = [];

pts_per_slice = 200;
pts_per_slice = 100;

spacing = 2.6;  %put vectors every this many microns

Lx = [-75 70];
Ly = [-65 65];
Lz = [-65 65];

% x = linspace(Lx(1),Lx(2),pts_per_slice);
% y = linspace(Ly(1),Ly(2),pts_per_slice);
% z = linspace(Lz(1),Lz(2),pts_per_slice);

x = Lx(1):spacing:Lx(2);
y = Ly(1):spacing:Ly(2);
z = Lz(1):spacing:Lz(2);

[X,Y,Z] = meshgrid(x,y,z); 
Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)]; 


%%

[body_tail] = combine_Meshes(Mesh, [1 2]); %%leave out transverse, wingtip since it's 2D

% [body_tail] = combine_Meshes(Mesh, [1 ]);

tic
[temp] = inside_mesh([Xlist(:) Ylist(:) Zlist(:)], body_tail);
disp(['is_inside took ',num2str(toc)]);
is_inside = NaN(size(Xlist));
is_inside(:) = temp;  %preserve dimensions of is_inside to match X, Y, Z

points = [Xlist(~is_inside) Ylist(~is_inside) Zlist(~is_inside)];


tic;
[field_vel] = field_velocity_mexed(points, Solutions.f{phase_ind}, Mesh, field_input);
disp(['field_velocity_mexed took ',num2str(toc)]);

U = NaN(size(X));  V = U;  W = U;
parfor i = 1:numel(X)
   
    point = [X(i) Y(i) Z(i)];
    ind = find(  points(:,1) == point(1) & points(:,2) == point(2) & points(:,3) == point(3) );
    U(i) = field_vel(ind,1);  V(i) = field_vel(ind,2);  W(i) = field_vel(ind,3);
end
%%


IN = inpolyhedron(body_tail.elems(:,1:3),body_tail.verts, x, y, z);

[X,Y,Z] = meshgrid(x,y,z);
inside.x = X(IN);  inside.y = Y(IN);  inside.z = Z(IN);
outside.x = X(~IN);  outside.y = Y(~IN);  outside.z = Z(~IN);

outside.x = X(~IN & X >= 0 & X <= 20);  outside.y = Y(~IN & X >= 0 & X <= 20);  outside.z = Z(~IN & X >= 0 & X <= 20);






nth = 1;
X2 = X(1:nth:end,1:nth:end,1:nth:end);
Y2 = Y(1:nth:end,1:nth:end,1:nth:end);
Z2 = Z(1:nth:end,1:nth:end,1:nth:end);
IN2 = IN(1:nth:end,1:nth:end,1:nth:end);
outside2.x = X2(~IN2);  outside2.y = Y2(~IN2);  outside2.z = Z2(~IN2);

outside2.x = X2(~IN2 & X2 >= 0 & X2 <= 20);  outside2.y = Y2(~IN2 & X2 >= 0 & X2 <= 20);  outside2.z = Z2(~IN2 & X2 >= 0 & X2 <= 20);

xlims = [9.2 10];
outside2.x = X2(~IN2 & X2 >= xlims(1) & X2 <= xlims(2));  outside2.y = Y2(~IN2 & X2 >= xlims(1) & X2 <= xlims(2));  outside2.z = Z2(~IN2 & X2 >= xlims(1) & X2 <= xlims(2));


% really need to remove inside velocities cause they can be very big....
U2 = U;  V2 = V;  W2 = W;
U2(IN) = NaN;  V2(IN) = NaN;  W2(IN) = NaN;

U0 = U;  V0 = V;  W0 = W;

U = U2;  V = V2;  W = W2;

% speed2 = sqrt(U2.^2+V2.^2+W2.^2);
%%
scale = 2;
close all
U = U2;  V = V2;  W = W2;

speed = sqrt(U2.^2+V2.^2+W2.^2);
speed_cutoff = quantile(speed(:),0.995);

too_big = speed > speed_cutoff;
U(too_big) = NaN;  V(too_big) = NaN;  W(too_big) = NaN;

% % vertical planes through transverse
% 
% x = near(6,X2); %bottom of groove  % 6
% lims = [x x; -Inf Inf; -Inf Inf];
% outside2.x = X2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );  
% outside2.y = Y2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% outside2.z = Z2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );  
% figure(1);  plot_mesh(Mesh); hold on;  axis tight;  light;
% cp = coneplot(X,Y,Z,U,V,W,outside2.x,outside2.y,outside2.z,scale);
% title('bottom of groove')
% 
% 
% 
x = near(9,X2); %middle of groove  % 9
lims = [x x; -Inf Inf; -Inf Inf];
outside2.x = X2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
outside2.y = Y2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
outside2.z = Z2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
figure(2);  plot_mesh(Mesh); hold on;  axis tight;  light;
cp = coneplot(X,Y,Z,U,V,W,outside2.x,outside2.y,outside2.z,scale,'quiver');
set(cp,'Color',[0 0 0]);
% hold on
% coord = [-20 0 -40];  vel = [100 0 0];  n = 3;
% xr = [coord(1)-1 coord(1) coord(1)+1]; yr = [coord(2)-1 coord(2) coord(2)+1]; zr = [coord(3)-1 coord(3) coord(3)+1];
% [Xr,Yr,Zr] = meshgrid(xr,yr,zr);
% 
% cpr = coneplot(Xr,Yr,Zr,repmat(vel(1),n,n,n),repmat(vel(2),n,n,n),repmat(vel(3),n,n,n),coord(1),coord(2),coord(3),scale,'quiver');
% cp = coneplot(X,Y,Z,U,V,W,scale,'nointerp');
title('middle of groove')
% 
% x = near(12.5,X2); %top of groove % 12.5
% lims = [x x; -Inf Inf; -Inf Inf];
% outside2.x = X2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% outside2.y = Y2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% outside2.z = Z2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% figure(3);  plot_mesh(Mesh); hold on;  axis tight;  light;
% cp = coneplot(X,Y,Z,U,V,W,outside2.x,outside2.y,outside2.z,scale);
% title('top of groove')
% 
% horizontal planes parallel to tail beating plane

z = near(0,Z2);  %middle of body, below tail % 0
lims = [-Inf Inf; -Inf Inf; z z];
outside2.x = X2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
outside2.y = Y2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
outside2.z = Z2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
figure(4);  plot_mesh(Mesh); hold on;  axis tight;  light;
cp = coneplot(X,Y,Z,U,V,W,outside2.x,outside2.y,outside2.z,scale,'quiver');

set(cp,'Color',[0 0 0]);
title('middle of body horizontal')
% 
% 
% z = near(10,Z2);  %bottom of tail % 10
% lims = [-Inf Inf; -Inf Inf; z z];
% outside2.x = X2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% outside2.y = Y2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% outside2.z = Z2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% figure(5);  plot_mesh(Mesh); hold on;  axis tight;  light;
% cp = coneplot(X,Y,Z,U,V,W,outside2.x,outside2.y,outside2.z,scale);
% title('bottom of tail')
% 
% 
% z = near(15,Z2);  %middle of tail  % 15
%  lims = [-Inf Inf; -Inf Inf; z z];
% outside2.x = X2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% outside2.y = Y2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% outside2.z = Z2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% figure(6);  plot_mesh(Mesh); hold on;  axis tight;  light;
% cp = coneplot(X,Y,Z,U,V,W,outside2.x,outside2.y,outside2.z,scale);
% title('middle of tail')
% 
% 
% 
% 
% z = near(20,Z2);  %top of tail % 20
%   lims = [-Inf Inf; -Inf Inf; z z];
% outside2.x = X2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% outside2.y = Y2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% outside2.z = Z2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
% figure(7);  plot_mesh(Mesh); hold on;  axis tight;  light;
% cp = coneplot(X,Y,Z,U,V,W,outside2.x,outside2.y,outside2.z,scale);
% title('top of tail')
% 
% 
%% 
% vertical, along tail axis, down through body
scale = 1/20;
U = U2;  V = V2;  W = W2;
speed = sqrt(U2.^2+V2.^2+W2.^2);
cutoff = 0.9925;  cutoff = 1;
speed_cutoff = quantile(speed(:),cutoff);
too_big = speed > speed_cutoff;
U(too_big) = NaN;  V(too_big) = NaN;  W(too_big) = NaN;
U = U*scale;  V = V*scale ;  W = W*scale;

y = near(0,Y2);  %through center of body, vertically, through tail % 0
  lims = [-Inf Inf; y y;   -Inf Inf;];
outside2.x = X2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
outside2.y = Y2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
outside2.z = Z2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
figure(989);  plot_mesh(Mesh); hold on;  axis tight;  light;
 cp = coneplot(X,Y,Z,U,V,W,outside2.x,outside2.y,outside2.z,0,'quiver');

set(cp,'Color',[0 0 0]);
% cp = coneplot(X,Y,Z,U,V,W,scale,'nointerp');
title('vertical through tail (scale bar = 100 \mum s^{-1})')
set(gca,'view',[0 0]);
hold on
coord = [-40 0 -40];  vel = [100*scale 0 0];  n = 3;
xr = [coord(1)-1 coord(1) coord(1)+1]; yr = [coord(2)-1 coord(2) coord(2)+1]; zr = [coord(3)-1 coord(3) coord(3)+1];
[Xr,Yr,Zr] = meshgrid(xr,yr,zr);

% scale arrow
% cpr = coneplot(Xr,Yr,Zr,repmat(vel(1),n,n,n),repmat(vel(2),n,n,n),repmat(vel(3),n,n,n),coord(1),coord(2),coord(3),0,'quiver');



%get the data from regular quiver
Ua = cp.UData;
Va = cp.VData;
Wa = cp.WData;
Xa = cp.XData;
Ya = cp.YData;
Za = cp.ZData;
delete(cp);

   arrowscale = 3  /  4;
start = [Xa Ya Za];
stop = start + arrowscale*[Ua Va Wa];
axis(axis)


[h] = arrow(start(~isnan(Va),:), stop(~isnan(Va),:),'length',6);

xlim([min(X2(:)) max(X2(:))]);  ylim([min(Y2(:)) max(Y2(:))]);  zlim([min(Z2(:)) max(Z2(:))]);



coord = [-40 0 -40];  vel = [100*scale 0 0];  n = 3;
xr = [coord(1)-1 coord(1) coord(1)+1]; yr = [coord(2)-1 coord(2) coord(2)+1]; zr = [coord(3)-1 coord(3) coord(3)+1];
[Xr,Yr,Zr] = meshgrid(xr,yr,zr);

% scale arrow
 cpr = coneplot(Xr,Yr,Zr,repmat(vel(1),n,n,n),repmat(vel(2),n,n,n),repmat(vel(3),n,n,n),coord(1),coord(2),coord(3),0,'quiver');
Uaa = cpr.UData;
Vaa = cpr.VData;
Waa = cpr.WData;
Xaa = cpr.XData;
Yaa = cpr.YData;
Zaa = cpr.ZData;
delete(cpr);


start = [Xaa Yaa Zaa];
stop = start + arrowscale*[Uaa Vaa Waa];
axis(axis)


[hr] = arrow(start(~isnan(Vaa),:), stop(~isnan(Vaa),:),'length',6);
% set(hr,'LineWidth',0.5);


xlim([min(X2(:)) max(X2(:))]);  ylim([min(Y2(:)) max(Y2(:))]);  zlim([min(Z2(:)) max(Z2(:))]);

%%
% middle of groove
scale = 0.15;
U = U2;  V = V2;  W = W2;
speed = sqrt(U2.^2+V2.^2+W2.^2);
cutoff = 0.9925;  cutoff = 1;
speed_cutoff = quantile(speed(:),cutoff);
too_big = speed > speed_cutoff;
U(too_big) = NaN;  V(too_big) = NaN;  W(too_big) = NaN;
U = U*scale;  V = V*scale ;  W = W*scale;

x = near(9,X2); %middle of groove  % 9
lims = [x x; -Inf Inf; -Inf Inf];
outside2.x = X2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
outside2.y = Y2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
outside2.z = Z2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
figure(2);  plot_mesh(Mesh); hold on;  axis tight;  light;
cp = coneplot(X,Y,Z,U,V,W,outside2.x,outside2.y,outside2.z,scale,'quiver');
set(cp,'Color',[0 0 0]);

% cp = coneplot(X,Y,Z,U,V,W,scale,'nointerp');
title('middle of groove (scale bar = 100 \mum s^{-1})')
set(gca,'view',[90 0]);
hold on



%get the data from regular quiver
Ua = cp.UData;
Va = cp.VData;
Wa = cp.WData;
Xa = cp.XData;
Ya = cp.YData;
Za = cp.ZData;
delete(cp);

   arrowscale = 3/4;
start = [Xa Ya Za];
stop = start + arrowscale*[Ua Va Wa];
axis(axis)


[h] = arrow(start(~isnan(Va),:), stop(~isnan(Va),:),'length',7);

xlim([min(X2(:)) max(X2(:))]);  ylim([min(Y2(:)) max(Y2(:))]);  zlim([min(Z2(:)) max(Z2(:))]);


coord = [x -40 -40];  vel = [0 100*scale 0];  n = 3;
xr = [coord(1)-1 coord(1) coord(1)+1]; yr = [coord(2)-1 coord(2) coord(2)+1]; zr = [coord(3)-1 coord(3) coord(3)+1];
[Xr,Yr,Zr] = meshgrid(xr,yr,zr);
% scale arrow
 cpr = coneplot(Xr,Yr,Zr,repmat(vel(1),n,n,n),repmat(vel(2),n,n,n),repmat(vel(3),n,n,n),coord(1),coord(2),coord(3),0,'quiver');

Uaa = cpr.UData;
Vaa = cpr.VData;
Waa = cpr.WData;
Xaa = cpr.XData;
Yaa = cpr.YData;
Zaa = cpr.ZData;
delete(cpr);


start = [Xaa Yaa Zaa];
stop = start + arrowscale*[Uaa Vaa Waa];
axis(axis)


[hr] = arrow(start(~isnan(Vaa),:), stop(~isnan(Vaa),:),'length',7);
% set(hr,'LineWidth',0.5);


xlim([min(X2(:)) max(X2(:))]);  ylim([min(Y2(:)) max(Y2(:))]);  zlim([min(Z2(:)) max(Z2(:))]);
%%
%%
% horizontal thru body
scale = 0.25;
U = U2;  V = V2;  W = W2;
speed = sqrt(U2.^2+V2.^2+W2.^2);
cutoff = 0.9925;  cutoff = 1;
speed_cutoff = quantile(speed(:),cutoff);
too_big = speed > speed_cutoff;
U(too_big) = NaN;  V(too_big) = NaN;  W(too_big) = NaN;
U = U*scale;  V = V*scale ;  W = W*scale;

z = near(0,Z2);  %middle of body, below tail % 0
lims = [-Inf Inf; -Inf Inf; z z];
outside2.x = X2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
outside2.y = Y2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
outside2.z = Z2(~IN2 & X2 >= lims(1,1) & X2 <= lims(1,2)  & Y2 >= lims(2,1) & Y2 <= lims(2,2)   & Z2 >= lims(3,1) & Z2 <= lims(3,2)  );
figure(4);  plot_mesh(Mesh); hold on;  axis tight;  light;
cp = coneplot(X,Y,Z,U,V,W,outside2.x,outside2.y,outside2.z,scale,'quiver');

set(cp,'Color',[0 0 0]);
title('middle of body horizontal (scale bar = 100 \mum s^{-1})')

set(gca,'view',[90 90]);
hold on



%get the data from regular quiver
Ua = cp.UData;
Va = cp.VData;
Wa = cp.WData;
Xa = cp.XData;
Ya = cp.YData;
Za = cp.ZData;
delete(cp);

   arrowscale = 3/4;
start = [Xa Ya Za];
stop = start + arrowscale*[Ua Va Wa];
axis(axis)


[h] = arrow(start(~isnan(Va),:), stop(~isnan(Va),:),'length',7);

xlim([min(X2(:)) max(X2(:))]);  ylim([min(Y2(:)) max(Y2(:))]);  zlim([min(Z2(:)) max(Z2(:))]);


coord = [59 -60 z];  vel = [0 100*scale 0];  n = 3;
xr = [coord(1)-1 coord(1) coord(1)+1]; yr = [coord(2)-1 coord(2) coord(2)+1]; zr = [coord(3)-1 coord(3) coord(3)+1];
[Xr,Yr,Zr] = meshgrid(xr,yr,zr);
% scale arrow
 cpr = coneplot(Xr,Yr,Zr,repmat(vel(1),n,n,n),repmat(vel(2),n,n,n),repmat(vel(3),n,n,n),coord(1),coord(2),coord(3),0,'quiver');

Uaa = cpr.UData;
Vaa = cpr.VData;
Waa = cpr.WData;
Xaa = cpr.XData;
Yaa = cpr.YData;
Zaa = cpr.ZData;
delete(cpr);


start = [Xaa Yaa Zaa];
stop = start + arrowscale*[Uaa Vaa Waa];
axis(axis)


[hr] = arrow(start(~isnan(Vaa),:), stop(~isnan(Vaa),:),'length',7);
% set(hr,'LineWidth',0.5);


xlim([min(X2(:)) max(X2(:))]);  ylim([min(Y2(:)) max(Y2(:))]);  zlim([min(Z2(:)) max(Z2(:))]);

%% subtract body velocity


r_cross_omega = cross( points - repmat(Mesh(1).refpoints(:,1)',size(points,1),1) ,  repmat(Answers.interaction_on.body_tail_transverse.kinematics.Omega',size(points,1),1)) ;
field_vel_bodyframe = field_vel - repmat(Answers.interaction_on.body_tail_transverse.kinematics.U',size(field_vel,1),1) - r_cross_omega;


figure;  [s,e] = plot_mesh(Mesh);  hold on;  set(gcf,'name','Vertical through transverse, body frame');
% middle of transverse
q = quiver3(  points(points(:,1)==8.75,1), points(points(:,1)==8.75,2),points(points(:,1)==8.75,3),field_vel_bodyframe(points(:,1)==8.75,1),field_vel_bodyframe(points(:,1)==8.75,2),field_vel_bodyframe(points(:,1)==8.75,3),6,'k');

figure;  [s,e] = plot_mesh(Mesh);  hold on;  set(gcf,'name','Horizontal through body center, below tail, parallel to tail beating plane, body frame');
% horizontal through middle of body, parallel to tail beating plane
q = quiver3(points( points(:,3)==0,1),points(points(:,3)==0,2),points(points(:,3)==0,3),field_vel_bodyframe(points(:,3)==0,1),field_vel_bodyframe(points(:,3)==0,2),field_vel_bodyframe(points(:,3)==0,3),10,'k');

figure;  [s,e] = plot_mesh(Mesh);  hold on;  set(gcf,'name','Horizontal through tail, parallel to tail beating plane, body frame');
% horizontal through tail, parallel to tail beating plane
q = quiver3(points( points(:,3)==15.25,1),points(points(:,3)==15.25,2),points(points(:,3)==15.25,3),field_vel_bodyframe(points(:,3)==15.25,1),field_vel_bodyframe(points(:,3)==15.25,2),field_vel_bodyframe(points(:,3)==15.25,3),10,'k');

figure;  [s,e] = plot_mesh(Mesh);  hold on; set(gcf,'name','Vertical along tail axis, down through body, body frame');
% vertical down through tail axis
q = quiver3(points( points(:,2)==0,1),points(points(:,2)==0,2),points(points(:,2)==0,3),field_vel_bodyframe(points(:,2)==0,1),field_vel_bodyframe(points(:,2)==0,2),field_vel_bodyframe(points(:,2)==0,3),20,'k');

%%






[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS

%X = 0;  Y = 0;  Z = 1;

[body_tail] = combine_Meshes(Mesh, [1 3]); %%leave out transverse since it's 2D
[temp] = inside_mesh([X(:) Y(:) Z(:)], body_tail);
disp('is_inside done');
is_inside = NaN(size(X));
is_inside(:) = temp;  %perserve dimensions of is_inside to match X, Y, Z

points = [X(~is_inside) Y(~is_inside) Z(~is_inside)];

tic;
[field_vel] = field_velocity_mexed(points, f, Mesh, field_input);
disp(['field_velocity_mexed took ',num2str(toc)]);

U = NaN(size(X));  V = U;  W = U;

U(~is_inside) = field_vel(:,1) * input.constants.multfactor / input.constants.mu;
V(~is_inside) = field_vel(:,2) * input.constants.multfactor / input.constants.mu;
W(~is_inside) = field_vel(:,3) * input.constants.multfactor / input.constants.mu;

Speed = sqrt(U.^2 + V.^2 + W.^2);
stopo


%%

Field.points = points;
Field.velocities = [U V W];
%%
x = linspace(-8,8,40); % x = 1.6;
y = linspace(-5,5,40);  
z = linspace(-5,5,40); 

% x = linspace(-2,3,19); % x = 1.6;
% y = linspace(-1.5,1.5,19);  
% z = linspace(-1.5,1.5,19); 
[X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
dist_tol = 0;
[idx]= rangesearch(Field.points,[X(:) Y(:) Z(:)],dist_tol,'distance','euclidean');
count = zeros(length(idx),1);
parfor i = 1:length(idx)
    count(i) = numel(idx{i});
end
points = [X(count == 0) Y(count == 0) Z(count == 0)];
[is_inside] = inside_mesh(points, body_tail);  %only compute at new points at least r distance from old
disp('is_inside done');

points = points(~is_inside,:); %further cut down size of points by excluding points inside mesh

tic;
[field_vel] = field_velocity_mexed(points, f, Mesh, field_input);
disp(['field_velocity_mexed took ',num2str(toc)]);
%%
%U = NaN(size(points));  V = U;  W = U;

U = field_vel(:,1);
V = field_vel(:,2) ;
W = field_vel(:,3);

Speed = sqrt(U.^2 + V.^2 + W.^2);

% concatenate new results with old
Field.points = [Field.points; points];
Field.velocities = [Field.velocities; [U V W] ];
Field.speed = sqrt(sum(Field.velocities.^2,2));

Field.interpolant = scatteredInterpolant(Field.points, Field.speed, 'natural', 'none');

bounds.x = [min(Field.points(:,1)) max(Field.points(:,1))];
bounds.y = [min(Field.points(:,2)) max(Field.points(:,2))];
bounds.z = [min(Field.points(:,3)) max(Field.points(:,3))];
npts = 300;
x = linspace(bounds.x(1),bounds.x(2),npts);
y = linspace(bounds.y(1),bounds.y(2),npts);
z = linspace(bounds.z(1),bounds.z(2),npts);

[X,Y,Z] = meshgrid(x,y,z);
Speed = Field.interpolant(X,Y,Z);

Field.interpolant.U = scatteredInterpolant(Field.points, Field.velocities(:,1), 'natural', 'none');
Field.interpolant.V = scatteredInterpolant(Field.points, Field.velocities(:,2), 'natural', 'none');
Field.interpolant.W = scatteredInterpolant(Field.points, Field.velocities(:,3), 'natural', 'none');

U = Field.interpolant.U(X,Y,Z);
V = Field.interpolant.V(X,Y,Z);
W = Field.interpolant.W(X,Y,Z);






%%
try, delete(sl); end;  try, delete(iso); end
xlim(bounds.x);  ylim(bounds.y);  zlim(bounds.z);
sl = slice(X,Y,Z,Speed,[ X(1,180,1)],[Y(151,1,1)],[ ],'linear');
set(sl,'facealpha',1);
set(sl,'edgealpha',0);
caxis([0.01 0.3]);
colormap jet
% fv = isosurface(X,Y,Z,Speed,0.12);
% iso = patch(fv);  set(iso,'facecolor','k','facealpha',0.5,'edgealpha',0.2);



%X(Y == Y(151,1,1) | X == X(1,180,1)),Y(151,1:inc:end,1:inc:end),Z(151,1:inc:end,1:inc:end),U(151,1:inc:end,1:inc:end),V(151,1:inc:end,1:inc:end),W(151,1:inc:end,1:inc:end)


%q = quiver3(X(151,1:inc:end,1:inc:end),Y(151,1:inc:end,1:inc:end),Z(151,1:inc:end,1:inc:end),U(151,1:inc:end,1:inc:end),V(151,1:inc:end,1:inc:end),W(151,1:inc:end,1:inc:end),1,'k','linewidth',1);

X_temp = X( X == X(1,169,1) | Z == Z(1,150,1));
Y_temp = Y( X == X(1,169,1) | X == X(1,150,1));
Z_temp = Z( X == X(1,169,1) | X == X(1,150,1));
U_temp = U( X == X(1,169,1) | X == X(1,150,1));
V_temp = V( X == X(1,169,1) | X == X(1,150,1));
W_temp = W( X == X(1,169,1) | X == X(1,150,1));

%W_temp(:) = 0;


X_temp = X(Y == Y(151,1,1) );
Y_temp = Y(Y == Y(151,1,1) );
Z_temp = Z(Y == Y(151,1,1) );
U_temp = (U(Y == Y(151,1,1) ));
V_temp = (V(Y == Y(151,1,1) ));
W_temp = (W(Y == Y(151,1,1) ));
V_temp(:) = 0;


%%
delete(q)
inc = 1;
q = quiver3(X_temp(1:inc:end),Y_temp(1:inc:end),Z_temp(1:inc:end),U_temp(1:inc:end),V_temp(1:inc:end),W_temp(1:inc:end),1.5,'k','linewidth',1);
figure(5)
%%
%q = quiver3(X(Y == Y(151,1,1) | X == X(1,180,1)),Y(Y == Y(151,1,1) | X == X(1,180,1)),Z(Y == Y(151,1,1) | X == X(1,180,1)),U(Y == Y(151,1,1) | X == X(1,180,1)),V(Y == Y(151,1,1) | X == X(1,180,1)),W(Y == Y(151,1,1) | X == X(1,180,1)),1,'k','linewidth',1);

lighting none
figure(5)