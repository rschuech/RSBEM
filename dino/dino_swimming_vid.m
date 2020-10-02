
% folder = 'C:\Users\rudi\Desktop\RD\dino_dumps_hairs_2_1.5_phases\';
% folder = 'C:\Users\rudi\Desktop\RD\hairs_2_1.5_newdumps\';
folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_coplanar_2_normal_1.5\';
base = 'Body-Transverse-Tail-Coplanar_Hairs-Normal_Top_Hairs-Normal_Bottom_Hairs';

folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail\';
base = 'Body-Transverse-Tail';

timestepping_file = [base,'_timestepping.mat'];
interp_file = [base,'_dump.mat'];

load([folder,interp_file],'Solutions');

% meshfolder =  'C:\Users\rudi\Desktop\RD\meshes_hairs_2_1.5_phases\';
meshfolder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\regular\final\';
parameters_file = [meshfolder, 'Parameters_step_1_time_0.000169836956522_phase_0.049087385212341.txt'];

do_align_path = true;
use_traction = false;

fat_tail = false;  coplanar_length = 2;  normal_length = 1.5;
videoname = [folder,base,'temp'];
savevid = true; % savevid = false;
do_hair_dots = false;
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

% names = {'Body' ,'Transverse','Tail'};


% dino_geom_parameters_huge_groove_uber_hairs;  % creates geom structure


% parameters_file = [folder_original , 'Parameters_step_',num2str(temp(1)), '_time_',num2str(temp(2),'%.15f'),'_phase_',num2str(temp(3),'%.15f'),'.txt'];

%     parameters_file = 'C:/Users/rudi/Desktop/RD/dino_mesh_automation/phases_hairs_4_1/original/Parameters_step_0_time_0.000000000000000_phase_0.000000000000000.txt';
% if ~exist(parameters_file,'file')
%     disp('parameters_file doesn''t exist; using file from first step');
%     parameters_file = [folder_original, 'Parameters_step_0_time_0.000000000000000_phase_0.000000000000000.txt'];
% end
dino_geom_parameters_from_file;


    
nthreads = 20;
%%
clear Inds Times Files Steps Phases dir
name = 'Metadata';
files = dir([meshfolder,name,'*']);
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

if use_traction
    nplaces = 10;
    Phases( ~ ismember( roundn(Phases,-nplaces) , roundn(Solutions.phase,-nplaces) )) = NaN;
end


%%


temp = load([folder,timestepping_file]);

if do_align_path
 [shift,angle,rotvec] = align_path(temp.fits,temp.timestepping_solution.y(1,:));
end


T = 1.5;  %max time for video  %30 for Body-Tail,    18 for Body-Tail-Transverse   65 for Body-Transverse-Coplanar_Hairs-Normal_Top_Hairs-Normal_Bottom_Hairs
  T = 4;
% full body revolution at around 6 sec or 1500 um

if T > max(temp.timestepping_solution.x)
    error('Problemo, requested T longer than length of timestepping solution');

end

[time, inds] = unique( temp.timestepping_solution.x( temp.timestepping_solution.x <= T) ) ;  %sometimes the ode45 output contains repeated time points
sol = temp.timestepping_solution.y( temp.timestepping_solution.x <= T , :) ;  sol = sol(inds,:);

interpolant = interp1(time, sol, 'spline' , 'pp');


dt = 1E-4;  % time step for video frames
dt = 0.025;
dt = 0.00025;
dt = 0.5;
 dt = 0.1;  % Body-Tail    Body-Tail-Transverse
  dt = 0.0005;  % Body-Transverse (refined enough to see direction of transverse wave)
%   dt = 0.05;  % Body-Transverse-Coplanar_Hairs-Normal_Top_Hairs-Normal_Bottom_Hairs
%    dt = .5;
dt = 0.001;

c = 0;  path_pts = [];
for t = 0:dt:T
    c = c+1;
    
    phase = t * geom.phase_speed;  %assumes initial angle = 0
    
    
    phase = mod( phase + 2*pi, 2*pi);  % shift periodic tail phase angle back to 0 - 2*pi
    
    [inds] = find( min(abs( phase - Phases)) == abs( phase - Phases) );
    if length(inds) > 1
        inds = inds(randperm(length(inds)));  
    end
    Step_closest = Steps(inds(1));
    Phase_closest = Phases(inds(1));
    metafile_closest = Files{inds(1)};
    
    % this seems needed for plotting traction, but has problems, main one
    % being when we haven't run the Fourier interpolation on  all mesh
    % files generated, and sometimes there are two equally closest phases
    % in Solutions and we pick the wrong one?
%     solutions_ind = find( min(abs( Phase_closest - Solutions.phase)) == abs( Phase_closest - Solutions.phase ) );
%     rand_inds = Solutions.rand_inds{solutions_ind};  % body tail transverse wingtip
    
    
    % it knows what Mesh to load and modify from metafile_closest
%     [Mesh, Vels, Parameters, Metadata] = interp_dino(t, meshfolder, metafile_closest, names, geom, nthreads, check_intersections,  rand_inds);
     [Mesh, Vels, Parameters, Metadata] = interp_dino(t, meshfolder, metafile_closest, names, geom, nthreads, check_intersections);
   
    
    if do_hair_dots
        u_in = linspace(0,Metadata.geom.transverse.u_max,150);
        h_in = [ Metadata.geom.transverse.h_max];
        
        [u_in, h_in] = ndgrid(u_in, h_in);  u_outer = u_in(:);  h_outer = h_in(:);
        
        
        u_in = linspace(0,Metadata.geom.transverse.u_max,800);
        h_in = [Metadata.geom.transverse.h_min];
        
        [u_in, h_in] = ndgrid(u_in, h_in);  u_inner = u_in(:);  h_inner = h_in(:);
        
        L_pts = length(u_inner);
        
        u_all = [u_inner; u_outer];  h_all = [h_inner; h_outer;];
        
        [edge_pts ] = transverse_hairs_parameterized(u_all, h_all, t, geom.transverse );
        
        Mesh(4).refpoints = [Mesh(4).refpoints  edge_pts];  %add edge points to wingtip refpoints
    end
    y = ppval(interpolant, t);
    
    
    scale = 1;
    for m = 1:length(Mesh)
        Mesh(m).verts = Mesh(m).verts * scale;
        Mesh(m).refpoints = Mesh(m).refpoints * scale;
    end
    
    Mesh2_temp = move_Mesh(Mesh,y);
    
    if do_align_path
    Mesh2 = shiftMesh(Mesh2_temp,shift);
    Mesh2 = rotateMesh(Mesh2,angle,rotvec);
    else
        Mesh2 = Mesh2_temp;
    end
    
    if use_traction
        Solutions_ind = find(roundn(Solutions.phase,-nplaces) == roundn(Phase_closest,-nplaces));   %this ind will be different from the ind going with metadata files because Solutions might have fewer phase evals and thus coarser resolution
        if numel(Solutions_ind) ~= 1
            problemo
        end
        
        f = Solutions.f{Solutions_ind};  %long column vector of solution tractions at global verts, alternating x y z components
        fmat = (reshape(f,3,[]))';  %row is x y z traction, col is global vert
        clear traction traction_mag
        for sf = 1:length(Mesh2)  %body transverse tail wingtip
            glob_inds = Mesh2(sf).indices.glob.vert;
            traction{sf} = fmat(glob_inds,:);
            traction_mag{sf} = sqrt(sum(traction{sf}.^2,2));
        end
    end
    
    
    path_pts(c,:) = Mesh2(1).Centroid;
    
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
      lims = [-15 240; -40 35; -40 35];
          
    hfig = figure(24);  hfig.Position = [ 415         185        1043         715]; % VideoWriter seems to not allow sizes much larger than this...
    clf
    
    
    
%     hold on
    if use_traction
        [s,e] = plot_mesh(Mesh2,[4 2 3 3 ], traction_mag, color_limits);
    else
        [s,e] = plot_mesh(Mesh2([1:3 ]),[2 2 1 2 2 2]);
    end
         set(e,'edgealpha',0.01);
    % set(s(4:6),'facealpha',0.75)
    
    hold on
    l = plot3(path_pts(:,1),path_pts(:,2),path_pts(:,3),'k--','linewidth',1);
    
    if do_hair_dots
        linel = plot3(Mesh2(4).refpoints(1,2:L_pts+1),Mesh2(4).refpoints(2,2:L_pts+1),Mesh2(4).refpoints(3,2:L_pts+1),'k-','linewidth',3);
        ptsl = plot3(Mesh2(4).refpoints(1,L_pts+2:end),Mesh2(4).refpoints(2,L_pts+2:end),Mesh2(4).refpoints(3,L_pts+2:end),'ko','markerfacecolor','k','markersize',4);
    end
    hold off
    
    xlim(lims(1,:));  ylim(lims(2,:)); zlim(lims(3,:));
  
    
    title(['Time = ',num2str(t,'%2.3f'),' s']);
    %     tl = text(-0, 20, -20,['Time = ',num2str(t),' s'],'fontsize',14);
    xlabel('\mum');
    ylabel('');
    zlabel('');
        set(gca,'YTick',[]);   set(gca,'ZTick',[]);
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
    set(s,'DiffuseStrength',0.7);
    
    
%     view([ -7.5742       19.035]);
%     view([-10.5 11.5]);
    view([-25 12]);
%     set(gca,'PlotBoxAspectRatio', [3.1167 1.3333 1]);
%     set(gca,'CameraTarget', [81.24 -10 -1.5655]);
%     set(gca,'CameraPosition', [-534.31 -857.24 164.3]);
%     set(gca,'CameraViewAngle', 4.76);
    %      end
    
%     set(s(3),'DiffuseStrength',1)
%     set(s(3),'SpecularStrength',1)
    
    
    %     set(gca,'CameraViewAngle',4.20880138025815,'DataAspectRatio',[1 1 1],...
    %     'PlotBoxAspectRatio',[1.07142857142857 1 1]);
    %
    
    
    % set(gca,'CameraPosition',...
    %     [-465.733272405367 -775.641925550281 351.124438779279],'CameraTarget',...
    %     [75.7550000421791 0.326530951544196 9.28245946211459],'CameraUpVector',...
    %     [0 0 1],'CameraViewAngle',3.80919250743569,'DataAspectRatio',[1 1 1],...
    %     'PlotBoxAspectRatio',[4.30798682874022 1 1.17349502586006],'YTick',[],...
    %     'ZTick',[]);
    
    
    % pause
    % axis off
    
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
%     pause
end

if savevid
    close(vidh);
end

