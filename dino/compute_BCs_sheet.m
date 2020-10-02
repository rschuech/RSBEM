folder = 'E:\Hull\dinoflagellate\meshes_stub_groups\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_stub_groups\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_grouped_fast\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail_nosphere2\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail_nosphere_tail_angle_0\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thick_tail_nosphere_tail_angle_0\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_bottom_flange\';
folder = 'E:\Hull\dinoflagellate\hairsheet3\';

do_vid = true;

if do_vid
    
    profile = 'MPEG-4';
    profile = 'Motion JPEG AVI';
    vidh = VideoWriter('E:\Hull\dinoflagellate\hairsheet3_vid',profile);
    
    vidh.FrameRate = 128/10; %50;
    vidh.Quality = 100; %1-100
    open(vidh);
    
end

input.performance.nthreads = 7;
store_constants = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear geom



geom.phase_speed = 2*pi*46;  %if lowest common period T = 4*pi s, covered in 2*pi rad, is 1/2 rad/sec for the entire beat cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% geom.sheet.lambda = 22.2;
geom.sheet.lambda = 4.7;
geom.sheet.nlambda = 3;
% geom.sheet.amp = [11.4 / 2  / 2];
geom.sheet.amp = [3.2  / 2];

geom.sheet.omega = 2 * pi * 46 ;   %rad/sec
% period = 2*pi rad / (2*pi*46  rad/sec) = 1/46 sec
% phase speed = 2*pi / period = 2*pi*46  rad/sec

geom.sheet.t_max = 2*pi * geom.sheet.nlambda;  %t value where hemispherical end cap is centered
% safety_factor = 1.5;  % t value at end of sphere must be less than safety_factor * geom.t_max
geom.sheet.hairlength = 0.3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Inds Files
names = {'Sheet'};
for n = 1:length(names)
    name = names{n};
    
    files = dir([folder,name,'*']);
    files = {files.name};
    %files = files(3:end);
    
    clear times
    for i = 1:length(files)
        ind1 = strfind(files{i},'time_') + 5;
        ind2 = strfind(files{i},'_phase') - 1;
        times(i) = str2double(files{i}(ind1:ind2));
    end
    
    [~,inds] = sort(times);
    Inds{n} = inds;
    Files{n} = files(inds);
    Times{n} = times(inds);
end


    Time = Times{1};



%%

for f = 1:length(Files{1})
    
    disp(['Step ',num2str(f),' out of ',num2str(length(Files{1}))]);
    
    time = Time(f);
    clear Mesh rand_inds
    
    [Mesh(1),rand_inds(1)] = load_mesh([folder,Files{1}{f}],[],[],'mesh');  %save rand vertex permutation for use later, so all future calculations use same exact mesh!
    
    
    Mesh(1).name = 'sheet';
    
    
    for si = 1:length(Mesh)
        Mesh(si).orientation = [1 0 0; 0 1 0; 0 0 1;]';
        Mesh(si).refpoints = [0 0 0]';
    end
    
    %tail indices might overlap body/transverse since they're not part of same
    %mesh
    %     Mesh(3).indices.orig.elem = Mesh(3).indices.orig.elem + max([Mesh(1).indices.orig.elem; Mesh(2).indices.orig.elem]);
    %     Mesh(3).indices.orig.vert = Mesh(3).indices.orig.vert + max([Mesh(1).indices.orig.vert; Mesh(2).indices.orig.vert]);
    %     Mesh(3).elems             = Mesh(3).elems             + max([Mesh(1).indices.orig.vert; Mesh(2).indices.orig.vert]);  %still using orig vert labels here, until renumber_Mesh below
    
    
    
    [Mesh] = global_inds(Mesh);
    [Mesh] = renumber_Mesh(Mesh);
    
    if store_constants
        
        temp_input.performance.nthreads = 7;
        
        temp_input.accuracy.integration_tol.area.abstol = 1E-9;
        temp_input.accuracy.integration_tol.area.reltol = 0.1;
        temp_input.accuracy.integration_tol.area.maxevals = Inf;
        
        temp_input.accuracy.integration_tol.centroid.abstol = 1E-12;
        temp_input.accuracy.integration_tol.centroid.reltol = 0.1;
        temp_input.accuracy.integration_tol.centroid.maxevals = Inf;
        
        temp_input.accuracy.integration_tol.volume.abstol = Inf;
        temp_input.accuracy.integration_tol.volume.reltol = Inf;  %sheet is 2D, volume meaningless
        temp_input.accuracy.integration_tol.volume.maxevals = Inf;
        
        [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
    end
    
    %     [is_intersected] = submesh_intersections(Mesh, geom.tail.radius / 8 ,[], true, 8);
    %     if is_intersected
    %         stopafra
    %     end
    
    %     names = {Mesh.name};
    %     ni = find(strcmp(names,'tail'));
    ni = 1;
    
    t_sol = NaN(1,Mesh(ni).n_vert);
    dist = t_sol;   sheet_vel = NaN(Mesh(ni).n_vert,3);  n_t = sheet_vel;
    pts = sheet_vel;  V_surfaces = pts;  V_us = pts; u_ns = pts;
    sheet_errors = NaN(length(t_sol),2);
    
    parfor(vert_i = 1:length(t_sol), 7)
      % for vert_i = 1:length(t_sol)
        %  vert_i / Mesh(ni).n_vert
        vert = Mesh(ni).verts(vert_i,:);
        %   vert = Surface_vis(ni).verts(vert_i,:);
        
        distance = @(t) sqrt(sum((vert - sheet_parameterized(t, time, geom.sheet)).^2));  %distance between vert and sheet centerline as function of t
        
%         [t_sol(vert_i), dist(vert_i), flag] = fminbnd(distance,0,geom.sheet.t_max*safety_factor,optimset('maxfunevals',1E15,'maxiter',1E15,'tolx',1E-14)); %limit upper bound by t_max so points on hemisphere will fall at this parameter value
        t_guess = vert(1) / geom.sheet.lambda * 2 * pi;  %solve for parameter t given x from vert, which shouldn't be far off actual t value
         [t_sol(vert_i), dist(vert_i), flag] = fminsearch(distance,t_guess,optimset('maxfunevals',1E15,'maxiter',1E15,'tolx',1E-10,'tolfun',1E-10)); %limit upper bound by t_max so points on hemisphere will fall at this parameter value
      
        if dist(vert_i) > geom.sheet.hairlength * 1.005  %allow 0.5% error
            stopafra
        end
        
       % dist(vert_i) = distance(t_sol(vert_i));
        
        [pt,der, u_n,pt_time(vert_i,:) ] = sheet_parameterized(t_sol(vert_i), time, geom.sheet);
        
%         if roundn(vert(1), -3) == 6.557
%             vert
%             pt
%             pause
%         end
        
        pts(vert_i,:) = pt;
        u_vert = vert - pt;  %vector from centerline pt to surface vert
        u_vert = u_vert / sqrt(sum(u_vert.^2));
        
        u_ns(vert_i,:) = u_n;
        
        u_verts(vert_i,:) = u_vert;
        
        %         dotprod = roundn( dot(u_vert,u_n/norm(u_n))  , -6);  %allow some numerical noise
%         if dist(vert_i) < 1E-3
%             dotprod = 0;
%         else
%             %           dotprod =  dot(u_vert,u_n/norm(u_n));
%             dotprod = roundn( dot(u_vert,u_n/norm(u_n))  , -4);  %allow some numerical noise
%         end
        
         dotprod =  dot(u_vert,u_n/norm(u_n));
         
        factor = sign(dotprod);   %fuxor cases where dotprod isn't 1 or -1 because of error in t_sol (i.e. when u_vert is 
        %so close to centerline that fminbnd can't figure out closest point
        %on centerline, and pt is actually correct.  In these cases, the
        %entire contribution of d(u_vert)/dt should be small anyway.
         
         
%         switch dotprod
%             case 1
%                 factor = 1;  % u_vert and u_n are in same direction, no problemo
%             case -1
%                 factor = -1;  %u_vert and u_n are in opposite directions, need to flip direction of u_vert_time, which is same as flipping analytical u_n (to make it actually align with u_vert)
%                 %this carries through and ends up just being multiplied by
%                 %u_vert_time
%             case 0
%                 factor = 0;
%             otherwise %numerical noise is too big
%                 stopafra
%                 
%         end
        
        %time derivative of normalized u_n multiplied by dist.  dist is constant for rigid hairs making up the sheet but need to normalize analytically
        u_vert_time(vert_i,:) = factor*dist(vert_i)*[-2*geom.sheet.amp*sin(-t_sol(vert_i)+geom.sheet.omega*time)*geom.sheet.omega*pi*geom.sheet.lambda^2/(4*geom.sheet.amp^2*cos(-t_sol(vert_i)+geom.sheet.omega*time)^2*pi^2+geom.sheet.lambda^2)^(3/2), ...
            -4*geom.sheet.lambda*geom.sheet.amp^2*cos(-t_sol(vert_i)+geom.sheet.omega*time)*sin(-t_sol(vert_i)+geom.sheet.omega*time)*geom.sheet.omega/(pi*(4*geom.sheet.amp^2*cos(-t_sol(vert_i)+geom.sheet.omega*time)^2+geom.sheet.lambda^2/pi^2)^(3/2)), ...
            0];
        
        % velocity of vert is vel of center line pt + vel of vector from centerline
        % pt to vert
        
        sheet_vel(vert_i,:) = pt_time(vert_i,:) + u_vert_time(vert_i,:);
        
    end
    factor  = 0.99;  %to make sure we don't look at points on flat end surface at t = tmax
    
    %     tail_error = abs([min(dist(t_sol < geom.tail.t_max*factor)) max(dist(t_sol < geom.tail.t_max*factor))] - geom.tail.radius);  %an estimate of error = difference between distances and actual radius
    %
    %       tail_errors(f,:) = tail_error;
    %
    %     if max(tail_error) > 0.5
    %         disp('stopafra')
    %         pause
    %     end
    
    
    
    
    %%
    %     t  = linspace(0,2*pi*geom.tail.nlambda*1.01,300);
    %     inc = 12;
    %     for Time = linspace(0,2,4000)
    %         [pt, vel, der] = tail_parameterized(t, Time, geom.tail);
    %
    %         figure(453)
    %         plot3(pt(:,1),pt(:,2),pt(:,3),'-','markerfacecolor','k');
    %         hold on
    %         %quiver3(pt(1:inc:end,1),pt(1:inc:end,2),pt(1:inc:end,3),vel(1:inc:end,1),vel(1:inc:end,2),vel(1:inc:end,3), 0.1);
    %         % q = quiver3(pt(1:inc:end,1),pt(1:inc:end,2),pt(1:inc:end,3),der(1:inc:end,1),der(1:inc:end,2),der(1:inc:end,3), 1);
    %         % set(q,'showarrowhead','off');
    %         [pt, vel, der] = tail_parameterized(geom.tail.t_max, Time, geom.tail);
    %         plot3(pt(:,1),pt(:,2),pt(:,3),'o','markerfacecolor','k');
    %
    %         hold off
    %         title(num2str(Time));
    %         grid on
    %         %ylim([-0.5 0.5]);
    %         view([0 90]);
    %        % pause
    %
    %         drawnow
    %
    %     end
    %%
    %     clear Pts
    %     tvec = linspace(0,geom.tail.t_max*1.05,300);
    %     for i = 1:length(tvec)
    %      [pt, vel, der] = tail_parameterized(tvec(i), 0, geom.tail);
    %      Pts(i,:) = pt;
    %     end
    %     hold on
    %     centre = plot3(Pts(:,1),Pts(:,2),Pts(:,3),'k-');
    %%
    
    % p = plot3(pts(:,1),pts(:,2),pts(:,3),'k.','markerfacecolor','k');
    %
    % inc = 80;
    % q = quiver3(pts(1:inc:end,1),pts(1:inc:end,2),pts(1:inc:end,3),V_surfaces(1:inc:end,1),V_surfaces(1:inc:end,2),V_surfaces(1:inc:end,3), 0.3);
    %
    % q = quiver3(pts(1:inc:end,1),pts(1:inc:end,2),pts(1:inc:end,3),V_us(1:inc:end,1),V_us(1:inc:end,2),V_us(1:inc:end,3), 0.3);
    
    
    %set(q,'showarrowhead','off');
    
    
    
    
    %% generate theoretical points more or less evenly spaced over transverse
    
    %     u = linspace(0,geom.transverse.u_max,5000);
    %     v = linspace(0,geom.transverse.v_max,100);
    %     c = 0;
    %     pts = NaN(length(u)*length(v),3);  U = NaN(length(u)*length(v),1); V = U;
    %     for i = 1:length(u)
    %         i/length(u)
    %
    %         for j = 1:length(v)
    %             c = c+1;
    %             U(c) = u(i);  V(c) = v(j);
    %
    %             pts(c,:) = transverse_parameterized(u(i),v(j), time, geom.transverse)';
    %         end
    %     end
    %
    %     hold on
    %     ph = plot3(pts(:,1),pts(:,2),pts(:,3),'k.');
    
    %% show closest points on theoretical surface to each vertex of actual mesh
    %     u = u_sol;  v = v_sol;  %u and v had better be same size since they match with verts
    %
    %     %c = 0;
    %     pts_fit = NaN(length(u),3);  U = NaN(length(u),1); V = U;
    %     parfor i = 1:length(u)
    %         %  i/length(u)
    %
    %         % c = c+1;
    %         %             U(c) = u(i);  V(c) = v(j);
    %
    %         pts_fit(i,:) = transverse_parameterized(u(i),v(i), time, geom.transverse)';
    %
    %     end
    %
    %     hold on
    %     ph = plot3(pts_fit(:,1),pts_fit(:,2),pts_fit(:,3),'k.');
    
    %%
    %     clear th
    %     for i = 1:size(pts_fit,1)
    %         th(i) = text(pts_fit(i,1),pts_fit(i,2),pts_fit(i,3),num2str(v_sol(i)),'fontsize',12);
    %     end
    %%
    
    %%
    clear Metadata
    Metadata.geom = geom;
    Metadata.time = time;
    Metadata.sheet.BCs = sheet_vel;
    Metadata.sheet.indices.orig.vert = Mesh(1).indices.orig.vert;
    Metadata.sheet.rand_inds = rand_inds(1).mesh.rand_inds;
    
    
    save([folder,'Metadata_',Files{1}{f}(6:end-4),'.mat'],'Metadata');
    
    
    
    
    
    figure(482)
    clf
    [s,e] = plot_mesh(Mesh, [1]);
    hold on
    
    
    %%
    
    try, delete(q1), end;   try, delete(q2), end;
    inc = 1;  ni = 1;
    q1 = quiver3(Mesh(ni).verts(1:inc:end,1),Mesh(ni).verts(1:inc:end,2),Mesh(ni).verts(1:inc:end,3),sheet_vel(1:inc:end,1),sheet_vel(1:inc:end,2),sheet_vel(1:inc:end,3), 1);
  
%     q1 = quiver3(Mesh(ni).verts(1:inc:end,1),Mesh(ni).verts(1:inc:end,2),Mesh(ni).verts(1:inc:end,3),u_ns(1:inc:end,1),u_ns(1:inc:end,2),u_ns(1:inc:end,3), 1);
    
    
    %q = quiver3(Surface_vis(ni).verts(1:inc:end,1),Surface_vis(ni).verts(1:inc:end,2),Surface_vis(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3));
    set(q1,'color','k');
     set(q1,'lineWidth',0.001)
    set(q1,'maxHeadSize',0.005)
    
    
%      xlim([-5       70]);  ylim([   -10     10]);  
 xlim([-1       15]);  ylim([   -2.5     2.5]);  
    
    % set(gca,'view',[  -24.321       47.549]);
    set(gcf,'position',[  582          31        1292         974]);
    % set(gca,'cameraViewAngle',4.0468);
    %     set(gca,'cameraPosition',[-128.15      -267.48       327.41]);
    %     set(gca,'cameraTarget',[-5.8668       3.0753       2.8308]);
    
    hold off
    
    title(['Step ',num2str(f),'       Time ',num2str(time)]);
    
    drawnow
    
   
    if do_vid
        
        hfig = gcf;
        orig_mode = get(hfig, 'PaperPositionMode');
        
        set(hfig, 'PaperPositionMode', 'auto');
        set(gcf,'position',[  582          31        1292         974]);
        cdata = hardcopy(hfig, '-Dopengl', '-r0');
        
        % Restore figure to original state
        
        set(hfig, 'PaperPositionMode', orig_mode); % end
        
        %  cdatas(:,:,:,n) = cdata;
        F = im2frame(cdata);
        
        writeVideo(vidh,F);
    end
    
    
end

if do_vid
    close(vidh)
end