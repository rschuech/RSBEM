
load('C:\Users\rudi\Desktop\RD\cleaner_dumps\curved_rod_AR1_1_AR2_0_tail_radius_0.03101752454497_amp_0.4953673272830541_lambda_3.472653357555359_nlambda_1.30595323359245_motorBC_torque_dump.mat');
Mesh0 = Mesh;

%%
videoname = 'C:\Users\rudi\Desktop\RD\bacteria_rotating_tail';
savevid = true;


try, close(vidh); end

if savevid
    profilee = 'Motion JPEG AVI';
    %   profile = 'Archival';
    vidh = VideoWriter(videoname,profilee);
    
    vidh.FrameRate = 30; %50;
    vidh.Quality = 100; %1-100
    open(vidh);
end

hfig = figure;

nthreads = 7;

c = 0; 
nrevs = 5;
degrees = repmat(linspace(0,360,120),1,nrevs);
for deg = degrees
    c = c+1;
    c/length(degrees)
    
 
   
     [Mesh(2)] = rotateMesh(Mesh0(2), [   0  0 deg*pi/180]' );  %rotate tail around x axis
        
      figure(hfig)
      cla
    [s,e] = plot_mesh(Mesh,[2 1]); light
    view([  -32         13.2]);
    axis off

    
  
    
   
    drawnow
    
    
    
    
    if savevid
        % Get CDATA from hardcopy using zbuffer
        
        % Need to have PaperPositionMode be auto
        
        orig_mode = get(hfig, 'PaperPositionMode');
        
        set(hfig, 'PaperPositionMode', 'auto');
        
        cdata = hardcopy(hfig, '-Dopengl', '-r0');
        
        % Restore figure to original state
        
        set(hfig, 'PaperPositionMode', orig_mode); % end
        
        %For the "OpenGL" renderer you can write a similar code. This technique will not work for the "painters" renderer.
        
        %Next, replace the use of GETFRAME from your code with IM2FRAME as follows:
        
        F = im2frame(cdata);
        
        writeVideo(vidh,F);
    end
    
end

if savevid
    close(vidh);
end

