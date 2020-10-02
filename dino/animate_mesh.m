

warning('off','MATLAB:hardcopy:DeprecatedHardcopyFunction');


% folder = 'E:\Hull\dinoflagellate\meshes_parallel2_3_all4\';
% folder = 'C:\Users\rudi\Desktop\RD\meshes_hairs_3.5_1.5\';
% folder = 'C:\Users\rudi\Desktop\RD\phases_hairs_3.5_1.5\';
folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\regular\cleaned\';

% videoname = 'C:\Users\rudi\Desktop\RD\BCs_3.5_1.5_transverse_tail_zoom_fixed';
videoname = 'C:\Users\rudi\Desktop\RD\pape main results\BCs_2_1.5_panels_testing';
savevid = true;

try, close(vidh); end

if savevid
    profile = 'MPEG-4';
    %   profile = 'Archival';
    vidh = VideoWriter(videoname,profile);
    
    vidh.FrameRate = 12; %50;
    vidh.Quality = 100; %1-100
    open(vidh);
end
%%

[Mesh_files] = get_dino_mesh_files(folder);

n_phase_pts = 128;
phase_list = linspace(0, 2*pi, n_phase_pts + 1);   % *always* 2 pi rad in a cycle
phase_list = phase_list(1:end-1);  % leave off last value which is same as first

for f = 1:length(phase_list)
    f
    tic
    phase = phase_list(f);
    
    temp_input.performance.nthreads = 20;
    temp_input.accuracy.integration_tol.area.abstol = 1E-6;
    temp_input.accuracy.integration_tol.area.reltol = 1000;
    temp_input.accuracy.integration_tol.area.maxevals = Inf;
    
    temp_input.accuracy.integration_tol.centroid.abstol = 1E-6;
    temp_input.accuracy.integration_tol.centroid.reltol = 1000;
    temp_input.accuracy.integration_tol.centroid.maxevals = Inf;
    
    temp_input.accuracy.integration_tol.volume.abstol = 1E-6;
    temp_input.accuracy.integration_tol.volume.reltol = 1000;
    temp_input.accuracy.integration_tol.volume.maxevals = Inf;
    
    temp_input.paths.datfolder = folder;
    temp_input.bugtype = 'dino';
    temp_input.potatohead = [1 1 1 1 1 1];
    temp_input.problemtype = 'freeswim';
    
    [Mesh, Metadata, matrix_props] = load_dino_mesh(phase, Mesh_files, temp_input);
    
    
%     [ Mesh , rotmat] = rotateMesh(Mesh,[0 -pi/2 0]);
     [ Mesh , rotmat] = rotateMesh(Mesh,[pi/2 0 0]);
    
    
    clear BCs
    for n = 1:length(Mesh)
        name = Mesh(n).name;
        [~,inds] = ismember(Mesh(n).indices.orig.vert,  Metadata.(name).indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs in Metadata to match current order in Mesh
        BCs.freeswim.(name) = Metadata.(name).BCs(inds,:);
        BCs.freeswim.(name) = ( rotmat * BCs.freeswim.(name)' )';
    end
    
    
    
    hfig = figure(582);
    clf
    
    subinds = {[1 5 6] , [1 5 6 2] , [1 5 6 3 4], [1:6]};
    
    for p = 1:4
        subtightplot(2,2,p,[0.02 0.02]);
    %     [s,e] = plot_mesh(Mesh(1:6),[2 3 2 2 1 2]);
    
    %     set(s(1),'FaceColor',[0 0.95 0]); %body
    %     set(s(6),'FaceColor',[0.1 0.1 0.8]); %transverse
    %     set(s(2),'FaceColor',[1 0.35 0.35],'facealpha',1); %coplanar
    %     set(s(3:4),'FaceColor',[0.75 0.75 1],'facealpha',1); %normal
    %     set(s(5),'FaceColor',[0.5 0.5 0.5]); %tail
    refines  = [3  4 3 3 2 3];
    [s,e] = plot_mesh(Mesh(subinds{p}),refines(subinds{p}));
%     set(s(1),'FaceColor','g'); %body
%     set(s(6),'FaceColor','r'); %transverse
%     set(s(5),'FaceColor','r'); %tail
%             set(s(2),'FaceColor','b','facealpha',1); %coplanar
%     set(s(3:4),'FaceColor','c','facealpha',1); %coplanar
    
    hold on
    transverse_ind = find(strcmp({Mesh.name},'Transverse'));
    tail_ind = find(strcmp({Mesh.name},'Tail'));
    coplanar_ind = find(strcmp({Mesh.name},'Coplanar_Hairs'));
    normal_top_ind = find(strcmp({Mesh.name},'Normal_Top_Hairs'));
    normal_bot_ind = find(strcmp({Mesh.name},'Normal_Bottom_Hairs'));
    
    hold off
    
%     try, delete(ph); end
    hold on
    u_in = linspace(0,Metadata.geom.transverse.u_max,150);
    
    cases = {    {'Transverse'} ,
        {'Coplanar_Hairs'} ,
        {'Normal_Top_Hairs' , 'Normal_Bottom_Hairs'} ,
        {'Coplanar_Hairs','Normal_Top_Hairs','Normal_Bottom_Hairs'}   };
    %     cases = {'Coplanar_Hairs'};
    for c = 1:length(cases{p})
        switch cases{p}{c}
            case 'Transverse'
                v_in = [Metadata.geom.transverse.v_max ];
            case 'Coplanar_Hairs'
                %                     h_in = [Metadata.geom.hairs.Coplanar_Hairs.h_min   Metadata.geom.hairs.(cases{c}).h_max];
                h_in = [Metadata.geom.hairs.Coplanar_Hairs.h_max ];
            case 'Normal_Top_Hairs'
                h_in = [Metadata.geom.hairs.(cases{p}{c}).h_min];
            case 'Normal_Bottom_Hairs'
                h_in = [Metadata.geom.hairs.(cases{p}{c}).h_max];
        end
        
        switch cases{p}{c}
            case 'Transverse'
                [u_in,v_in] = ndgrid(u_in,v_in);  u_in = u_in(:);  v_in = v_in(:);
                [pt] = transverse_parameterized(u_in, v_in, Metadata.time, Metadata.geom.transverse );
            otherwise
                [u_in,h_in] = ndgrid(u_in,h_in);  u_in = u_in(:);  h_in = h_in(:);
                [pt] = transverse_hairs_parameterized(u_in, h_in, Metadata.time, Metadata.geom.transverse , cases{p}{c});
        end
        pt = rotmat * pt;
        
        plot3(pt(1,:),pt(2,:),pt(3,:),'ko','markerfacecolor','k','markersize',2);
    end
    hold off
    
    set(e,'edgeAlpha',0)
    %     set(s,'facealpha',1)
    

    %zoomed in to groove and transverse
%     xlim([   -12       12]);
%     ylim([   20        40]);
%     zlim([   -Inf       Inf]);
    %showing entire cell
    xlim([ -18.1877   18.1874]);
    ylim([ -11.2695   47.6002]);
    zlim([ -18.4496   20.4210]);
    view([64.8434   26.5250]);
    
    
    
%     l = light('position',[-1 0 -1]);
    
     l = light('position',[1 -1 1]);
     
     
    set(s,'DiffuseStrength',0.7);
    %     tex = text(20,15,-20,file,'interpreter','none','fontsize',15);
    str = ['Step ',num2str(f),'     Time ',num2str(Metadata.time),'     Phase ',num2str(phase)];
    %      tex = text(-19,5,17,str,'interpreter','none','fontsize',15);
    annotation('textbox','Position',[0.46 0.96 0.12 0.036],'String',['step ',num2str(f-1)],'LineStyle','none');
    %     set(gcf,'position',[427          65        1301         913]);
    
%     set(gcf,'position',[680         107        1180         871]);
%     set(gcf,'position',[-1583         169        1129         947]);
    set(gcf,'position',[748          45        1129         947]);
    
    % zoomed in
%     set(gca,'CameraPosition',...
%         [-510.969312765543 -134.903208869965 13.8172335915585],'CameraTarget',...
%         [-3.11237710715562 2.29984741219147 31.0387862366723],'CameraUpVector',...
%         [0 0 1],'CameraViewAngle',1.5793,'DataAspectRatio',[1 1 1],...
%         'PlotBoxAspectRatio',[1.26792342254347 1 2.0000244188318]);
    
    %zoomed out
    % set(gca,'CameraPosition',...
    %     [-507.917685337482 -137.167971025864 9.24245183443854],'CameraTarget',...
    %     [-0.0607496790947704 0.0350852562922495 26.4640044795523],'CameraUpVector',...
    %     [0 0 1],'CameraViewAngle',2.36876254790011,'DataAspectRatio',[1 1 1],...
    %     'PlotBoxAspectRatio',[1.26792342254347 1 2.0000244188318]);
    
    %zoomed out showing tail beat
    % set(gca,'CameraPosition',...
    %     [-330.663718222777 -407.457608594561 -4.21928545291752],'CameraTarget',...
    %     [2.1369993699174 -0.342352840995511 19.0353598356348],'CameraUpVector',...
    %     [0 0 1],'CameraViewAngle',5.32663520672285,'DataAspectRatio',[1 1 1],...
    %     'PlotBoxAspectRatio',[1.26792342254347 1 2.0000244188318]);
    
    %     xlim([ -27.709       27.709]);
    %     ylim([   -20.449       20.449]);
    %     zlim([   -44.11       37.686]);
    axis off
    
    end
    
    f / length(phase_list)
    
    drawnow
    
%     return
    
    
    if savevid
 
        cdata = print(hfig, '-RGBImage','-r150'); %r0 works but is allegedly only 96 DPI.  150 works but anything higher seems to break.
      
        F = im2frame(cdata);
        
        writeVideo(vidh,F);
    end
   
    toc
end

if savevid
    close(vidh);
end


