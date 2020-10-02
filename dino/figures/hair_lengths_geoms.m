
files = {'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_coplanar_3_normal_1.5\Body-Transverse-Tail-Coplanar_Hairs-Normal_Top_Hairs-Normal_Bottom_Hairs_dump.mat',...
    'C:\Users\rudi\Desktop\RD\pape main results\body_tail_transverse_coplanar_3\Body-Transverse-Tail-Coplanar_Hairs_dump.mat',...
    'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_normal_1.5\Body-Transverse-Tail-Normal_Top_Hairs-Normal_Bottom_Hairs_dump.mat',...
    'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail\Body-Transverse-Tail_dump.mat'};

names = {'3 1.5', '3 0', '0 1.5', '0 0'};


for f = 1:length(files)
    file = files{f}; name = names{f};
    
    load(file,'Mesh');
    
    
    figure(94);
    clf;
    set(gcf,'Position',[ 680         113        1223         985]);
    [s,e] = plot_mesh(Mesh,2);
    set(e,'EdgeAlpha',0.1);
    L = light; L.Position = [0 -1 0];
    view([-25 20]);
    xlim([-11.5 48]); ylim([-21.5 21.5]); zlim([-21.5 21.5]);
    % grid off
    axis off
    
    export_fig(gcf,['C:\Users\rudi\Desktop\RD\pape main results\figures\hair geoms ',name,'.png'],'-r400','-transparent');
    
end

