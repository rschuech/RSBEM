clear sphere_eff sphere_temporal curved_eff curved_temporal

load('C:\Hull\methods_fig_dumps\optimized tails comparison\curved_rod_AR1_1_AR2_0_tail_radius_0.03101752454497_amp_0.461328125_lambda_3.296875_nlambda_1.3296875_motorBC_torque_dump.mat','Mesh','Metadata');
sphere_eff.Mesh = Mesh; sphere_eff.Metadata = Metadata;
load('C:\Hull\methods_fig_dumps\optimized tails comparison\curved_rod_AR1_1_AR2_0_tail_radius_0.03101752454497_amp_0.2592667655289663_lambda_3.908794398745886_nlambda_3.548452451793616_motorBC_torque_dump.mat','Mesh','Metadata');
sphere_temporal.Mesh = Mesh; sphere_temporal.Metadata = Metadata;
load('C:\Hull\methods_fig_dumps\optimized tails comparison\curved_rod_AR1_10_AR2_0.5_tail_radius_0.03101752454497_amp_0.836328125_lambda_5.8046875_nlambda_1.3765625_motorBC_torque_dump.mat','Mesh','Metadata');
curved_eff.Mesh = Mesh;  curved_eff.Metadata = Metadata;
load('C:\Hull\methods_fig_dumps\optimized tails comparison\curved_rod_AR1_10_AR2_0.5_tail_radius_0.03101752454497_amp_0.6090091278597636_lambda_8.050746774978951_nlambda_1.690616806304947_motorBC_torque_dump.mat','Mesh','Metadata');
curved_temporal.Mesh = Mesh;  curved_temporal.Metadata = Metadata;


figure(368);  clf; set(gcf,'Position',[    -1882          30        1786         922]);
[s1,e1] = plot_mesh(shiftMesh(sphere_eff.Mesh,[sphere_eff.Metadata(1).geom.radius 0 0]),[4 2]);
hold on
[s2,e2] = plot_mesh(shiftMesh(sphere_temporal.Mesh(2),[sphere_temporal.Metadata(1).geom.radius 0 0]),2);

[s3,e3] = plot_mesh(shiftMesh(rotateMesh(curved_eff.Mesh,[0 0 0]),[curved_eff.Metadata(1).geom.radius -1.5 0]) ,[4 2]);
hold on
[s4,e4] = plot_mesh(shiftMesh(rotateMesh(curved_temporal.Mesh(2),[0 0 0]),[curved_temporal.Metadata(1).geom.radius -1.5 0]),2 );
hold off
axis tight;  axis off;  light;
%%
e1(1).EdgeAlpha = 0.3;  % edges body sphere
e3(1).EdgeAlpha = 0.3;  % edges body curved rod
e1(2).EdgeAlpha = 0.0;  e2(1).EdgeAlpha = 0.0; e3(2).EdgeAlpha = 0.0;  e4(1).EdgeAlpha = 0.0;
s1(2).FaceColor = 'r';  s2(1).FaceColor = 'b'; s3(2).FaceColor = 'r';  s4(1).FaceColor = 'b';


% exportfig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\optimized tails comparison.png','format','png','color','rgb','resolution',800,'bounds','loose')