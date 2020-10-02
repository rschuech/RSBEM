debug_mode = false;  % plot just one bug to see details of path
align_paths = true;
savevid = false;
clear Meshes Interpolants Fits Solutions


if debug_mode
    videoname = ['C:\Users\rudi\Desktop\RD\',interp_dumps(1:end-24)];
else
    videoname = 'C:\Users\rudi\Desktop\RD\bacteria_race9';
end


T = 0.6;

folder = 'C:\Users\rudi\Desktop\RD\swimming dumps\';
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
    max(time)
    sol = sol(time <= T,:); refpoints = refpoints(time <= T,:);   time = time(time <= T);
    
    interpolant = interp1(time, sol, 'spline' , 'pp');
    interpolant2 = interp1(time, refpoints, 'spline' , 'pp');
    
    Meshes{d} = dump.Mesh;
    Interpolants{d} = interpolant;
    Interpolants2{d} = interpolant2;
    Fits(d) = dump.fits;
    Solutions(d) = dump.timestepping_solution;
    
end

Meshes = Meshes([4 1 3 5 2 6]);
Interpolants = Interpolants([4 1 3 5 2 6]);
Interpolants2 = Interpolants2([4 1 3 5 2 6]);
Fits = Fits([4 1 3 5 2 6]);
Solutions = Solutions([4 1 3 5 2 6]);
interp_dumps = interp_dumps([4 1 3 5 2 6])';

if debug_mode
    dump_ind = 4; % which bug to do
    
    Meshes = Meshes(dump_ind); Interpolants = Interpolants(dump_ind); Interpolants2 = Interpolants2(dump_ind); Fits = Fits(dump_ind); Solutions = Solutions(dump_ind); interp_dumps = interp_dumps{dump_ind};
end

clear shift angle rotvec
for d = 1:length(Meshes)
    [shift(d,:),angle(d,:),rotvec(d,:)] = align_path(Fits(d),Solutions(d).refpoint(1,:));
end
%%
if debug_mode
    lims = [-Inf Inf; -Inf Inf; -Inf Inf;];
    shifts = 0;
else
    lims = [-8 17.5; -5 5; -24 3;];
    shifts = linspace(0,-20,6);
end

dt = 0.001;

% dt = 0.005;





try, close(vidh); end

if savevid
    profilee = 'Motion JPEG AVI';
    %   profile = 'Archival';
    vidh = VideoWriter(videoname,profilee);
    
    vidh.FrameRate = 30; %50;
    vidh.Quality = 100; %1-100
    open(vidh);
end






c = 0;  path_pts = [];  tail_angle = [];  swimming_axis = []; tail_axis = [];
for t = 0:dt:T
    c = c+1;
    
    
    hfig = figure(23);
    clf
    
    
    
    aligned_paths = [];
    for d = 1:length(Meshes)
    
        n = size(Solutions(d).refpoint,1);
        aligned_path = Solutions(d).refpoint(:,1:3)' + repmat(shift(d,:)',1,n);
        rotmat = rotate_arbitrary_vector( rotvec(d,:), angle(d,:));
        aligned_path = rotmat * aligned_path;
        aligned_path(3,:) = aligned_path(3,:)  + shifts(d);
        aligned_paths{d} = aligned_path;
        
        plot3(aligned_path(1,:),aligned_path(2,:),aligned_path(3,:),'r-');
        hold on
        
    end
    
    for d = 1:length(Meshes)
      
        interpolant = Interpolants{d};
        Mesh = Meshes{d};
        
        y = ppval(interpolant, t);
        
        [Mesh(2)] = rotateMesh(Mesh(2), [0 0 y(7)]' );  %rotate tail around x axis
        %         y(3) = y(3) + shifts(d);
        
        Mesh = move_Mesh(Mesh,y);
        
        if align_paths
            Mesh = shiftMesh(Mesh,shift(d,:));
            Mesh = rotateMesh(Mesh,angle(d,:),rotvec(d,:));
        end
        
        Mesh = shiftMesh(Mesh,[0 0 shifts(d)]);
        
        
        path_pts{d}(c,:) = Mesh(1).refpoints(:,1);
        %         tail_angle{d}(c) = acos(dot( Mesh(1).orientation(:,1)  ,  [1 0 0]'  ));
        
        %         tail_axis{d}(c,:) = Mesh(1).orientation(:,1);
        
        % rotvec_temp = cross(Mesh(1).orientation(:,1),[1 0 0]');  rotvec_temp = rotvec_temp / sqrt(sum(rotvec_temp.^2));
        
        %   rotmat = rotate_arbitrary_vector( rotvec_temp,  tail_angle{d}(c));
        
        rotmat = Mesh(1).orientation';  % the coordinate transformation matrix from fixed frame to body frame is simply the transpose of the body frame basis vectors
        
        swimming_axis{d}(c,:) = rotmat * [1 0 0]';  swimming_axis{d}(c,:) = swimming_axis{d}(c,:)/ sqrt(sum(swimming_axis{d}(c,:).^2));
        
        %   if d == 5 && c == 300
        %       pause
        %   end
        %         [s,e] = plot_mesh(Mesh,[2 1]);
        %            [s,e] = plot_mesh(Mesh,[1 0]);
        %         %     set(e,'edgealpha',0);
        %         set(e,'edgealpha',0.1)
        %     set(s,'ambientStrength',0.7)
        %          light
        
        hold on
        l = plot3(path_pts{d}(:,1),path_pts{d}(:,2),path_pts{d}(:,3),'k--','linewidth',1);
        axis equal
        hold on
    end
    
    
    hold off
    
    xlim(lims(1,:));  ylim(lims(2,:)); zlim(lims(3,:));
    title(['Time = ',num2str(t),' s']);
    %     tl = text(-0, 20, -20,['Time = ',num2str(t),' s'],'fontsize',14);
    xlabel('\mum');
    ylabel('');
    zlabel('');
    set(gca,'YTick',[]);   set(gca,'ZTick',[]);
    %     grid off
    %      if t == 0
    %     lig = light('position',[1 -1 0]);
    
    view([0 0]);
    
    drawnow
    %     sdfsdfsdf
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






