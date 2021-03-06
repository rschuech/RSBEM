
folders = {'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_coplanar_2_normal_1.5\'};

bases = {'Body-Transverse-Tail-Coplanar_Hairs-Normal_Top_Hairs-Normal_Bottom_Hairs'};


nthreads = 20;

% time step for video frames

   dt = 0.25;
% dt = 0.001;  % for single dino
dt = 0.002; % race #1 2 or  body tail

T = 1.5;  %max time for video  %30 for Body-Tail,    18 for Body-Tail-Transverse   65 for Body-Transverse-Coplanar_Hairs-Normal_Top_Hairs-Normal_Bottom_Hairs
T = 4;
T = Inf; % race 1 or 2
T = temp.timestepping_solution.x(end);
% T = 12; % body tail
% full body revolution at around 6 sec or 1500 um

% meshfolder =  'C:\Users\rudi\Desktop\RD\meshes_hairs_2_1.5_phases\';
meshfolders = {'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\regular\final\',...
     'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\stopped transverse\'};

 % appears to be some problem with Metadata for stopped_transverse folder, 
 % because the plotting doesn't work, but the mesh files seem fine according to quick mesh plot.
 % Metadata files seem to have zero velocity BCs for transverse etc so perhaps indices are fuxored
 % try redoing compute BCs?

parameters_file = [meshfolder, 'Parameters_step_1_time_0.000169836956522_phase_0.049087385212341.txt'];


fat_tail = false;  coplanar_length = 2;  normal_length = 1.5;

videoname = 'C:\Users\rudi\Desktop\RD\pape main results\videos\turning2';
savevid = true; %savevid = false;

check_intersections = false;  %don't bother if we're already sure there are no intersection problems

try, close(vidh); end

if savevid
    vidprofilee = 'MPEG-4';
    %   vidprofile = 'Archival';
    vidh = VideoWriter(videoname,vidprofilee);
    
    vidh.FrameRate = 130; %50;  % 130 seems to be the max allowed?
    vidh.Quality = 75; %1-100
    open(vidh);
end


% names = {'Body',  'Tail','Transverse', 'Coplanar_Hairs','Normal_Top_Hairs','Normal_Bottom_Hairs'};
names = {'Body' ,'Transverse','Tail', 'Coplanar_Hairs','Normal_Top_Hairs','Normal_Bottom_Hairs'};


%%
clear Inds Times Files Steps Phases dir
name = 'Metadata';
files = dir([meshfolders{1},name,'*']);
files = {files.name};

clear times  phases
for i = 1:length(files)
    ind0 = strfind(files{i},'step_') + 5;
    ind01 = strfind(files{i},'_time') - 1;
    ind1 = strfind(files{i},'time_') + 5;
    ind2 = strfind(files{i},'_phase') - 1;
    ind3 = strfind(files{i},'phase_') + 6;
    ind4 = strfind(files{i},'.mat') - 1;
    
    times(i) = str2double(files{i}(ind1:ind2));
    steps(i) = str2double(files{i}(ind0:ind01));
    phases(i) = str2double(files{i}(ind3:ind4));
end

[~,inds] = sort(times);

Files = files(inds);
Times = times(inds);
Steps = steps(inds);
Phases = phases(inds);

%%
clear geoms interpolants shifts angles rotvecs
for f = 1:length(folders)
    folder = folders{f};  base = bases{f};
    
%     timestepping_file = [base,'_timestepping.mat'];
%     interp_file = [base,'_dump.mat'];
    
%     load([folder,interp_file],'Solutions');
    dino_geom_parameters_from_file;  geoms{f} = geom;
    
%     temp = load([folder,timestepping_file]);
    
    [time, inds] = unique( temp.timestepping_solution.x( temp.timestepping_solution.x <= T) ) ;  %sometimes the ode45 output contains repeated time points
    sol = temp.timestepping_solution.y( temp.timestepping_solution.x <= T , :) ;  sol = sol(inds,:);
    
    interpolants{f} = interp1(time, sol, 'spline' , 'pp');
end

hfig = figure(24); 

c = 0;  path_pts = []; clear leading_x
for t = 0:dt:T
    c = c+1;
    
    % VideoWriter seems to not allow sizes much larger than this...
%        hfig = figure(24);  hfig.Position = [ 415         185        1043         715]; %single dino
  
 hfig.Position = [    573         530        1300         460]; %race 1 or 2
%  hfig.Position = [658   233   945   651]; % body tail
hfig.Position = [605         530        1268         460]; % body centered tail
hfig.Position = [1052         196         821         794];
               clf
    
    geom = geoms{1};
    %%
    phase = t * geom.phase_speed;  %assumes initial angle = 0
%    phase
    
    phase = mod( phase + 2*pi, 2*pi);  % shift periodic tail phase angle back to 0 - 2*pi
    
    [inds] = find( min(abs( phase - Phases)) == abs( phase - Phases) );
    if length(inds) > 1
        inds = inds(randperm(length(inds)));
    end
    Step_closest = Steps(inds(1));
    Phase_closest = Phases(inds(1));
    metafile_closest = Files{inds(1)};
    
    % it knows what Mesh to load and modify from metafile_closest
    %    [Mesh, Vels, Parameters, Metadata] = interp_dino(t, meshfolder, metafile_closest, names, geom, nthreads, check_intersections,  rand_inds);
   if t*input.phase_speed >= piecewise_interpolant.phase_breaks(1) && t*input.phase_speed <= piecewise_interpolant.phase_breaks(2)
       meshfolder = meshfolders{1};
   elseif t*input.phase_speed >= piecewise_interpolant.phase_breaks(2) && t*input.phase_speed <= piecewise_interpolant.phase_breaks(3)
       meshfolder = meshfolders{2};
   elseif t*input.phase_speed >= piecewise_interpolant.phase_breaks(3)
       meshfolder = meshfolders{1};
   end
       
      meshfolder = meshfolders{1};
    
    [Mesh, Vels, Parameters, Metadata] = interp_dino(t, meshfolder, metafile_closest, names, geom, nthreads, check_intersections);
    
    if c == 1
        Mesh_0 = Mesh;
    end
    
    %%
    % this will be the "master mesh" with all components
%     Mesh(3) = shiftMesh(Mesh(3),geom.tail.shift);
%     n_pts = 15;  theta = linspace(0,2*pi,n_pts);  rho = repmat(body_radius,1,n_pts);    % body_radius should be defined from dino_geom_parameters_from_file.m
%     [y_dots,z_dots] = pol2cart(theta,rho); x_dots = repmat(40,1,n_pts);
%     Mesh(1).refpoints = [Mesh(1).refpoints   [x_dots; y_dots; z_dots;]  ];
    
%     Mesh(1).refpoints = [Mesh(1).refpoints  [42.25 0 15]'  [50 0 15]'  ];
    
    scale = 1;
    for m = 1:length(Mesh)
        Mesh(m).verts = Mesh(m).verts * scale;
        Mesh(m).refpoints = Mesh(m).refpoints * scale;
    end
    
    Mesh0 = Mesh;
    Mesh_00 = Mesh_0;
    
    for f = 1:length(folders)
        Mesh = Mesh0;
%         Mesh(1).Centroid
        y = ppval(interpolants{f}, t);
        
        removes = [];
        for mm = 1:length(Mesh)
            if ~contains(bases{f}, Mesh(mm).name)
                removes(end+1) = mm;
            end
        end
        Mesh(removes) = [];
        
        Mesh2_temp = move_Mesh(Mesh,y);
        Mesh_00 = move_Mesh(Mesh_00,y);
  
        Mesh2 = Mesh2_temp;
  
        path_pts(c,:,f) = Mesh2(1).Centroid;
        
        %     l = 250;
        lims = [-35 3000; -40 40; -30 150];  %zoomed out for helical path
        %     lims = [-40 90; -50 50; -80 60] ;
        lims = [-20 1000; -40 40;  -60 60]; % Body Tail Transverse
        lims = [5 400; -20 20;  -20 20]; %  % Body-Transverse
        lims = [-500 5500; -700 700;  -700 700];
        %         lims = [-23.4 166.1; -20.5 23.6; -20.8 30.9];  %zoomed out for helical path
        lims = [-10 60; -80 -0; 300 380];
        
        lims = [-12 175; -50 30; -35 25]; % whole dino?
        lims = [-Inf Inf; -Inf Inf; -Inf Inf];
        lims = [-12 110; -60 30; -65 25];
        lims = [-15 240; -40 35; -40 35]; % single dino e.g. body-transverse-tail
        lims = [-60 900; -125 150; -180 70]; % race # 1
        lims = [-60 900; -125 150; -200 95]; % race # 2
%           lims = [-60 70; -125 150; -75 50]; %tail only
         lims = [-25 150; -125 150; -150 30]; % body centered tail
        
        
        hold on
        
        
          if t*input.phase_speed >= piecewise_interpolant.phase_breaks(1) && t*input.phase_speed <= piecewise_interpolant.phase_breaks(2)
       [s,e] = plot_mesh(Mesh2([1 2 3 4 5 6]),[2 2 1 2 2 2]); % race 1 or 2 or  body tail
   elseif t*input.phase_speed >= piecewise_interpolant.phase_breaks(2) && t*input.phase_speed <= piecewise_interpolant.phase_breaks(3)
      [s,e] = plot_mesh(Mesh2(3),1); % race 1 or 2 or  body tail
      hold on
         [s2,e2] = plot_mesh(Mesh_00([1 2 4 5 6]),[2 2 2 2 2 2]); % race 1 or 2 or  body tail
        s = [s s2];  e = [e e2];
   elseif t*input.phase_speed >= piecewise_interpolant.phase_breaks(3)
      [s,e] = plot_mesh(Mesh2([1 2 3 4 5 6]),[2 2 1 2 2 2]); % r
   end
        
        
        
     
        set(e,'edgealpha',0.01);
    set(s,'DiffuseStrength',0.7);


        % set(s(4:6),'facealpha',0.75)
        
        hold on
%         dots = plot3(Mesh2(1).refpoints(1,2:end),Mesh2(1).refpoints(2,2:end),Mesh2(1).refpoints(3,2:end),'k.','MarkerSize',4);
%        arr = arrow3(Mesh2(1).refpoints(:,2)',Mesh2(1).refpoints(:,3)' , 'k-1');
        hold off
%         pause
    end
    
    hold on
    for pp = 1:size(path_pts,3)
        l(pp) = plot3(path_pts(:,1,pp),path_pts(:,2,pp),path_pts(:,3,pp),'k--','linewidth',1);
        
    end
    hold off
    
        xlim(lims(1,:));  ylim(lims(2,:)); zlim(lims(3,:));
    
    
    title(['Time = ',num2str(t,'%2.3f'),' s']);
    %     tl = text(-0, 20, -20,['Time = ',num2str(t),' s'],'fontsize',14);
    xlabel('\mum');
    ylabel('\mum');
    zlabel('\mum');
%     set(gca,'YTick',[]);   set(gca,'ZTick',[]);
    grid off
    %      if t == 0
    lig = light('position',[1 -1 0]);
    %     view([-76 2.4]);
    %     view([ -33.806       26.638]);
    
    %     view([-90 0]);
    %     view([-62 6.4]);  % Body-Tail
    %     view([-86.4 2.4]); % Body Tail Transverse
    %      view([-36 9]);% Body-Transverse
    %       view([  90 0 ]);

    
    
    %     view([ -7.5742       19.035]);
    %     view([-10.5 11.5]);
    view([-25 12]); % for individual dino, e.g. body-transverse-tail
    view([0 0]);  % race
   
    

    
    % pause
    % axis off
%     set(gca,'Position',[0.04 0.02 0.96 0.96]);
    drawnow
    
    
    if savevid
        % Get CDATA from hardcopy using zbuffer
        
        % Need to have PaperPositionMode be auto
        
        %         orig_mode = get(hfig, 'PaperPositionMode');
        
        %         set(hfig, 'PaperPositionMode', 'auto');
        
        %         cdata = hardcopy(hfig, '-Dopengl', '-r0');
        cdata = print(hfig, '-RGBImage');
        % Restore figure to original state
        
        %         set(hfig, 'PaperPositionMode', orig_mode); % end
        
        %For the "OpenGL" renderer you can write a similar code. This technique will not work for the "painters" renderer.
        
        %Next, replace the use of GETFRAME from your code with IM2FRAME as follows:
        
        F = im2frame(cdata);
        
        writeVideo(vidh,F);
    end
%         pause
end

if savevid
    close(vidh);
end

