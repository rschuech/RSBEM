
savevid = true;

type = 'rotation';

videoname = ['C:\Users\rudi\Desktop\RD\ellipsoid','_',type];

try, close(vidh); end

if savevid
    profilee = 'Motion JPEG AVI';
    %   profile = 'Archival';
    vidh = VideoWriter(videoname,profilee);
    
    vidh.FrameRate = 30; %50;
    vidh.Quality = 100; %1-100
    open(vidh);
end

boxlims = [-2.5 2.5; -2.5 2.5; -2.5 2.5];

shifts = linspace(0,3,100);

angles = linspace(0,4* 2*pi,4 * 150);
fact = 10;
linestyle = 'k:';

hfig = figure(1234);
set(gcf,'position',[    35         508        1846         596]);
p = panel();
p.pack(1,3);
p.margin = [0 0 0 0];
 p(1,1).select();
p(1,2).select();
p(1,3).select();

set(gcf,'color','w');
switch type
    case 'rotation'
looped = angles;
    case 'translation'
        looped = shifts;
end
for i = 1:length(looped)
    switch type
    case 'rotation'
 i/length(angles)
    case 'translation'
       i/length(shifts)
end
   
 %  subplot(1,3,1)
   p(1,1).select();
      cla
          switch type
    case 'rotation'
 Mesh2 = move_Mesh(Mesh,[[ 0 0 0] angles(i) 0 0 0]');
    case 'translation'
        Mesh2 = move_Mesh(Mesh,[[ shifts(i) 0 0] 0 0 0 0]');
end
   
   [s,e] = plot_mesh(Mesh2);  
set(e,'EdgeAlpha',0.1);
 

axis equal; 
  light
   lighting phong
    
    set(s,'ambientStrength',0.7)
        xlim([boxlims(1,1) boxlims(1,2)])
    ylim([boxlims(2,1) boxlims(2,2)])
    zlim([boxlims(3,1) boxlims(3,2)])
        set(gca,'view',[-45 10]);
    axis off
   % box on    
        hold on
    plot3(fact*[boxlims(1,:)],[0 0],[0 0],linestyle)
      plot3([0 0],fact*[boxlims(2,:)],[0 0],linestyle)
      plot3([0 0],[0 0],fact*boxlims(3,:),linestyle)
    hold off
    
     %  subplot(1,3,2)
        p(1,2).select();
         cla
                   switch type
    case 'rotation'
 Mesh2 = move_Mesh(Mesh,[[ 0 0 0] 0 angles(i) 0  0]');
    case 'translation'
        Mesh2 = move_Mesh(Mesh,[[ 0 shifts(i) 0] 0 0 0 0]');
end
   
     [s,e] = plot_mesh(Mesh2);  
set(e,'EdgeAlpha',0.1);
   

axis equal; 
  light
    lighting phong
    
    set(s,'ambientStrength',0.7)
    xlim([boxlims(1,1) boxlims(1,2)])
    ylim([boxlims(2,1) boxlims(2,2)])
    zlim([boxlims(3,1) boxlims(3,2)])
        set(gca,'view',[-45 10]);
    axis off
       % box on
            hold on
    plot3(fact*[boxlims(1,:)],[0 0],[0 0],linestyle)
      plot3([0 0],fact*[boxlims(2,:)],[0 0],linestyle)
      plot3([0 0],[0 0],fact*boxlims(3,:),linestyle)
    hold off
        
        
      % subplot(1,3,3)
        p(1,3).select();
         cla
                            switch type
    case 'rotation'
 Mesh2 = move_Mesh(Mesh,[[ 0 0 0] 0 0 angles(i)  0]');
    case 'translation'
        Mesh2 = move_Mesh(Mesh,[[ 0 0 shifts(i) ] 0 0 0 0]');
end
   
     [s,e] = plot_mesh(Mesh2);  
set(e,'EdgeAlpha',0.1);
   
axis equal; 
  light
    lighting phong
    
    set(s,'ambientStrength',0.7)
    xlim([boxlims(1,1) boxlims(1,2)])
    ylim([boxlims(2,1) boxlims(2,2)])
    zlim([boxlims(3,1) boxlims(3,2)])
        set(gca,'view',[-45 10]);
   axis off
    %box on
    hold on
    plot3(fact*[boxlims(1,:)],[0 0],[0 0],linestyle)
      plot3([0 0],fact*[boxlims(2,:)],[0 0],linestyle)
      plot3([0 0],[0 0],fact*boxlims(3,:),linestyle)
    hold off
    
    
drawnow

if savevid            
          orig_mode = get(hfig, 'PaperPositionMode');        
            set(hfig, 'PaperPositionMode', 'auto');         
            cdata = hardcopy(hfig, '-Dopengl', '-r0');     
            set(hfig, 'PaperPositionMode', orig_mode); % end
            F = im2frame(cdata);
            
            writeVideo(vidh,F);
        end

end

if savevid
    close(vidh);
end