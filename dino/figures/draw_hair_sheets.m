load('C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_coplanar_2_normal_1.5\Body-Transverse-Tail-Coplanar_Hairs-Normal_Top_Hairs-Normal_Bottom_Hairs_dump.mat');
geom = Metadata_temp.geom;  time = 0;

figure(11);  clf;  set(gcf,'Position',[977    49   931   947]);
%     pan = panel();
% pan.pack(3,2);
% pan.margin = repmat(0.01,1,4);

n_hairs = 204 * 8;  % Salome says length of transverse axoneme is 204 um and Miyasaka 2004 assumed 8 hairs / um

lims = [24.25 31.75; 2.5 14; 2.5 20];

submeshes = {6,[6 2],6,[6 3 4],6,[6 2 3 4]};
% panel_inds = [1 1; 1 2; 2 1; 2 2; 3 1; 3 2;];

clear hair_handles
for p = 1:6
    
    subtightplot(3,2,p,[0.02 0.02]);
    
    
    %  pan(panel_inds(p,1),panel_inds(p,2)).select();
    
    
    
    % 88.2634   72.7296
    [s,e] = plot_mesh(Mesh(submeshes{p}),4);  light;  %set(gca,'View',[146 33]);
    set(e,'EdgeAlpha',0.1);
    set(gca,'View',[90 90]);
    
    switch p
        case 1
            hair_handles{p} = sketch_hairs('Coplanar_Hairs',n_hairs,geom,time);
            %            [ph] = sketch_axoneme(geom,time);
        case 3
            hair_handles{p} = sketch_hairs('Normal_Top_Hairs',n_hairs,geom,time);  hair_handles{p} = sketch_hairs('Normal_Bottom_Hairs',n_hairs,geom,time);
            %        [ph] = sketch_axoneme(geom,time);
        case 5
            hair_handles{p} = sketch_hairs('Coplanar_Hairs',n_hairs,geom,time);
            hair_handles{p} = sketch_hairs('Normal_Top_Hairs',n_hairs,geom,time);  hair_handles{p} = sketch_hairs('Normal_Bottom_Hairs',n_hairs,geom,time);
            %    [ph] = sketch_axoneme(geom,time);
    end
    [ph] = sketch_axoneme(geom,time);
    axis tight
    xlim(lims(1,:)); ylim(lims(2,:)); zlim(lims(3,:));
    %     axis off
    set(gca,'XLabel',[]);  set(gca,'YLabel',[]);  set(gca,'ZLabel',[]);
    set(gca,'XTick',[]);  set(gca,'YTick',[]);  set(gca,'ZTick',[]);
end

% export_fig(gcf,['C:\Users\rudi\Desktop\RD\pape main results\figures\hair sheets actual density','.png'],'-r400','-transparent');



function [ph] = sketch_hairs(hair_type,n_hairs,geom,time)
u = linspace(0,geom.transverse.u_max,n_hairs);
v = linspace(geom.hairs.(hair_type).h_min,geom.hairs.(hair_type).h_max,2);
c = 0;
%         pts = NaN(length(u)*length(v),3);
clear pts
U = NaN(length(u)*length(v),1); V = U;
for i = 1:length(u)
    clear temp
    cc = 0;
    for j = 1:length(v)
        c = c+1;  cc = cc + 1;
        U(c) = u(i);  V(c) = v(j);
        %pts(c,:) = transverse_parameterized(u(i),v(j), time, geom.transverse)';
        %                 [pts(c,:)] = transverse_hairs_parameterized(u(i),v(j), time, geom.transverse , hair_type);
        temp(:,cc) = transverse_hairs_parameterized(u(i),v(j), time, geom.transverse , hair_type); % rows = x y z, cc = 1, 2
        
    end
    pts.x(:,i) = temp(1,:);  pts.y(:,i) = temp(2,:);  pts.z(:,i) = temp(3,:);
end

hold on
ph = plot3(pts.x,pts.y,pts.z,'k-','linewidth',1.5);
end



function [ph] = sketch_axoneme(geom,time)
u = linspace(0,geom.transverse.u_max,500);
v = 0;  hair_type = 'Coplanar_Hairs';

%         pts = NaN(length(u)*length(v),3);
clear pts

for i = 1:length(u)
    
    
    
    %pts(c,:) = transverse_parameterized(u(i),v(j), time, geom.transverse)';
    %                 [pts(c,:)] = transverse_hairs_parameterized(u(i),v(j), time, geom.transverse , hair_type);
    pts(:,i) = transverse_hairs_parameterized(u(i),v, time, geom.transverse , hair_type); % rows = x y z, cc = 1, 2
    
    
    %     pts.x(:,i) = temp(1,:);  pts.y(:,i) = temp(2,:);  pts.z(:,i) = temp(3,:);
end

hold on
ph = plot3(pts(1,:),pts(2,:),pts(3,:),'k-','linewidth',4);
end


