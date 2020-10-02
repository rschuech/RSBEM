 load('C:\Users\rudi\Desktop\RD\efficiency_dumps_pole2pole\curved_rod_AR1_5_AR2_0.5_tail_radius_0.03101752454497_amp_0.4762219362959685_lambda_3.405836583485502_nlambda_1.387773625462044_motorBC_torque_dump.mat')
%%
 figure(45); clf;  [s,e] = plot_mesh(Mesh,2);
set(gcf,'Color','w');

   [~,sphere_centers, centerline_center] = calc_pole_coords(Metadata);
   rear_to_front_dir = sphere_centers.head - sphere_centers.tail;  rear_to_front_dir = rear_to_front_dir / sqrt(sum(rear_to_front_dir.^2));
   
   line_pt = [sphere_centers.tail(1) sphere_centers.tail(2) 1]';
   
            pole2pole_line = [line_pt' rear_to_front_dir'];
           
            lh = drawLine3d(pole2pole_line); hold on;
         
            
            
%             set(s,'facealpha',1)
%               set(e,'edgealpha',0.1)
%               lig = light;
              
              light
lighting phong

set(s,'ambientStrength',0.5);
set(e,'edgealpha',0.075);
set(s,'facealpha',0.4);
set(s,'facecolor',repmat(0.7,1,3));
              
              
              lh.LineStyle = '--';
              lh.Color = 'k';
              lh.LineWidth = 1.5;
              
            sc(1) =   plot3(sphere_centers.head(1),sphere_centers.head(2),1,'ko','MarkerFaceColor','k','MarkerSize',5);
             sc(2) =   plot3(sphere_centers.tail(1),sphere_centers.tail(2),1,'ko','MarkerFaceColor','k','MarkerSize',5);
             
             ylim([-5 2.5]);
             xlim([-3 2]);
            
axis off

