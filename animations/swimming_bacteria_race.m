debug_mode = false;  % plot just one bug to see details of path
align_paths = true;
savevid = false;

movie_type = 'race';  % 'race' of several shapes or 'single' of one shape zoomed on path

clear Meshes Interpolants Fits Solutions


if debug_mode
    videoname = ['C:\Users\rudi\Desktop\RD\',interp_dumps(1:end-24)];
else
    videoname = 'C:\Users\rudi\Desktop\RD\bacteria_race_pole2pole3';
%      videoname = 'C:\Users\rudi\Desktop\RD\path_zoom';
end


% T = 0.6;
T = 1;
finish_line_x = 15;

folder = 'C:\Users\rudi\Desktop\RD\swimming dumps\pole2pole\';  clear dir
interp_dumps = dir([folder,'*dump.mat']);
interp_dumps = {interp_dumps.name};

timestepping_dumps = dir([folder,'*timestepping.mat']);
timestepping_dumps = {timestepping_dumps.name};


for d = 1:length(interp_dumps)
    interp_dump = load([folder,interp_dumps{d}]);
    timestepping_dump = load([folder,timestepping_dumps{d}]);
    dump = interp_dump;
    fields = fieldnames(timestepping_dump);
    for f = 1:length(fields)
        dump.(fields{f}) = timestepping_dump.(fields{f});
    end
    
    Avg_Power = dump.avg_omega * dump.input.tail.motor_torque;  %one of these should be constant and the other varying unless I eventually implement a constant power motor condition
    Adj_Speed = sqrt( dump.fits.converged.speed.^2 .* dump.input.constants.power  ./ Avg_Power);
    speedup_factor = Adj_Speed  /  dump.fits.converged.speed  % actual speeds / original speeds that applies to everything in the problem
    
    [time, inds] = unique( dump.timestepping_solution.x ) ;  %sometimes the ode113 output contains repeated time points
    sol = dump.timestepping_solution.y(inds,:);
    refpoints = dump.timestepping_solution.refpoint(inds,:);
    
    time = time / speedup_factor;  %time gets crunched by this factor - things take less time if the speeds are actually faster
%     max(time)
if max(time) < T
    error('Simulated time < requested T');
end

    sol = sol(time <= T,:); refpoints = refpoints(time <= T,:);   time = time(time <= T);
    
    interpolant = interp1(time, sol, 'spline' , 'pp');
    interpolant2 = interp1(time, refpoints, 'spline' , 'pp');
    
    Meshes{d} = dump.Mesh;
    Interpolants{d} = interpolant;
    Interpolants2{d} = interpolant2;
%     if d == 6
%         stoap
%     end
    temp = dump.fits;
    if isfield(temp,'avg_swimming_axis')
        temp = rmfield(temp,'avg_swimming_axis');
    end
    Fits(d) = temp;
    Solutions(d) = dump.timestepping_solution;
    
end

inds_order = [1 5 4 6 2 3 ];
% inds_order = [1 5 7 6 2 3 ];

Meshes = Meshes(inds_order);
Interpolants = Interpolants(inds_order);
Interpolants2 = Interpolants2(inds_order);
Fits = Fits(inds_order);
Solutions = Solutions(inds_order);
interp_dumps = interp_dumps(inds_order)';


if debug_mode
    dump_ind = 4; % which bug to do
    
    Meshes = Meshes(dump_ind); Interpolants = Interpolants(dump_ind); Interpolants2 = Interpolants2(dump_ind); Fits = Fits(dump_ind); Solutions = Solutions(dump_ind); interp_dumps = interp_dumps{dump_ind};
end

clear shift angle rotvec
for d = 1:length(Meshes)
    [shift(d,:),angle(d,:),rotvec(d,:)] = align_path(Fits(d),Solutions(d).refpoint(1,:));
end



% clear leading_x
% for i = 1:length(Meshes)
%     Mesh = Meshes{i};
% %      if align_paths
% %             Mesh = shiftMesh(Mesh,shift(i,:));
% %             Mesh = rotateMesh(Mesh,angle(i,:),rotvec(i,:));
% %         end
%     leading_x(i) = max(Mesh(1).verts(:,1));
% end
% race_shifts = max(leading_x) - leading_x;
% for i = 1:length(Meshes)
%     Meshes{i} = shiftMesh(Meshes{i},[race_shifts(i) 0 0]);
% end



%%
if debug_mode
    lims = [-Inf Inf; -Inf Inf; -Inf Inf;];
    shifts = 0;
else
      switch movie_type
        case 'race'
          %     lims = [-8 17.5; -5 5; -24 3;];
     lims = [-10.5 16; -5 5; -21.5 1.5];
    
        case 'single'
            lims = [-8 16; -2 2; -6 -2];
      end
    shifts = linspace(0,-20,6);

end

dt = 0.001;

%   dt = 0.025;

dt = 0.2;



try, close(vidh); end

if savevid
%     profilee = 'Motion JPEG AVI';
    profilee = 'MPEG-4';
    %   profile = 'Archival';
    vidh = VideoWriter(videoname,profilee);
    
    vidh.FrameRate = 30; %50;
    vidh.Quality = 95; %1-100
    resolution = '-r0';
%     resolution = '-r180';  % doesn't seem to look much different, but
%     doubles file size
    open(vidh);
end

race_align_leading_edges;  % computes race_shifts, the x shifts to apply to each bug to align the leading edges after all the other shifts, rotations are done

c = 0;  path_pts = [];  tail_angle = [];  swimming_axis = []; tail_axis = [];  
for t = 0:dt:T
    c = c+1;
    
    
    hfig = figure(27);  set(hfig,'Position',[  720         155        1028         900]);
    clf
    
     try, delete(leading_edges); end;
    
    aligned_paths = [];
    for d = 1:length(Meshes)
        
        n = size(Solutions(d).refpoint,1);
        aligned_path = Solutions(d).refpoint(:,1:3)' + repmat(shift(d,:)',1,n);
        rotmat = rotate_arbitrary_vector( rotvec(d,:), angle(d,:));
        aligned_path = rotmat * aligned_path;
        aligned_path(3,:) = aligned_path(3,:)  + shifts(d);
        aligned_paths{d} = aligned_path;
        
        %         plot3(aligned_path(1,:),aligned_path(2,:),aligned_path(3,:),'r-');
        %         hold on
        
    end
    
    clear leading_edges
    switch movie_type
        case 'race'
            D = 1:length(Meshes);
        case 'single'
            D = 2;
    end
    for d = 4 %D
        
        interpolant = Interpolants{d};
        Mesh = Meshes{d};
        
        y = ppval(interpolant, t);
    [t    y']
        
        %         [Mesh(2)] = rotateMesh(Mesh(2), [0 0 y(7)]' );  %rotate tail around x axis
        
        
        
        [Mesh(2)] = rotateTail(Mesh(2), y(7));
        
        
        
        Mesh = move_Mesh(Mesh,y);
        
        if align_paths
            Mesh = shiftMesh(Mesh,shift(d,:));
            Mesh = rotateMesh(Mesh,angle(d,:),rotvec(d,:));
        end
        
        Mesh = shiftMesh(Mesh,[0 0 shifts(d)]);
        
        
        Mesh = shiftMesh(Mesh,[race_shifts(d) 0 0]);
        
        
        path_pts{d}(c,:) = Mesh(1).refpoints(:,1);
        %         tail_angle{d}(c) = acos(dot( Mesh(1).orientation(:,1)  ,  [1 0 0]'  ));
        
        %         tail_axis{d}(c,:) = Mesh(1).orientation(:,1);
        
        % rotvec_temp = cross(Mesh(1).orientation(:,1),[1 0 0]');  rotvec_temp = rotvec_temp / sqrt(sum(rotvec_temp.^2));
        
        %   rotmat = rotate_arbitrary_vector( rotvec_temp,  tail_angle{d}(c));
        
        rotmat = Mesh(1).orientation';  % the coordinate transformation matrix from fixed frame to body frame is simply the transpose of the body frame basis vectors
        
%         swimming_axis{d}(c,:) = rotmat * [1 0 0]';  swimming_axis{d}(c,:) = swimming_axis{d}(c,:)/ sqrt(sum(swimming_axis{d}(c,:).^2));
        
         leading_edges(d) = max(Mesh(1).verts(:,1));
%            leading_x = max(Mesh(1).verts(:,1));
%           leading_edges(d) = plot3(repmat(leading_x,1,2),zeros(1,2),[-20 0],'k--');
          
        %   if d == 5 && c == 300
        %       pause
        %   end
        %         [s,e] = plot_mesh(Mesh,[2 1]);
        [s,e] = plot_mesh(Mesh,[1 0]);
        %         %     set(e,'edgealpha',0);
        set(e,'edgealpha',0.1)
        set(s,'ambientStrength',0.4)
        %                  light
        
        hold on
        l = plot3(path_pts{d}(:,1),path_pts{d}(:,2),path_pts{d}(:,3),'k-','linewidth',0.5);
        axis equal
        hold on
        
        
    end
    
    
    hold off

    xlim(lims(1,:));  ylim(lims(2,:)); zlim(lims(3,:));
    set(gca,'fontsize',12);
    title(['Time = ',num2str(t,'%.3f'),' s']);
    %     tl = text(-0, 20, -20,['Time = ',num2str(t),' s'],'fontsize',14);
            switch movie_type
        case 'race'
            xlabel('\mum','fontsize',14);
    ylabel('');
    zlabel('');
    set(gca,'YTick',[]);   set(gca,'ZTick',[]);
      view([0 0]);
        case 'single'
            view(gca,[-74.8 6.4]);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[6 1 1],'YTick',...
    [],'ZTick',[]);
           axis off
    end
   
    set(gca,'Position',[  0.0400    0.0673    0.9250    0.9077]);
    %     grid off
%     if t == 0
        lig = light('position',[1 -1 0]);
%           lig = light('position',[20 -10 -12]);
        lighting gouraud
%     end
  
    
    drawnow
    %     sdfsdfsdf
    if savevid
     
 
        orig_mode = get(hfig, 'PaperPositionMode');
        
        set(hfig, 'PaperPositionMode', 'auto');
        
%         cdata = hardcopy(hfig, '-Dopengl', '-r0');
        cdata = print('-RGBImage','-opengl',resolution);
        % Restore figure to original state
        
        set(hfig, 'PaperPositionMode', orig_mode); % end
        
        %For the "OpenGL" renderer you can write a similar code. This technique will not work for the "painters" renderer.
        
        %Next, replace the use of GETFRAME from your code with IM2FRAME as follows:
        
        F = im2frame(cdata);
        
        writeVideo(vidh,F);
    end
    
if any(leading_edges > finish_line_x)
    break
end
    
    
end

if savevid
    close(vidh);
end






