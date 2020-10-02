folder = 'E:\Hull\dinoflagellate\meshes_stub_groups\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_stub_groups\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_grouped_fast\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail_nosphere2\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail_nosphere_tail_angle_0\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thick_tail_nosphere_tail_angle_0\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_bottom_flange\';

% folder = 'E:\Hull\dinoflagellate\meshes_transverse_hairs_bot\';
% folder = 'E:\Hull\dinoflagellate\meshes_transverse_hairs_wider_parallel\';
folder = 'E:\Hull\dinoflagellate\meshes_transverse_hairs_parallel2_3_all\';
folder = 'E:\Hull\dinoflagellate\meshes_parallel2_3\';

outfolder = 'E:\Hull\dinoflagellate\meshes_parallel2_3_modded\';

input.performance.nthreads = 40;
folder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_3\';

folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing2\';

outfolder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_3_modded\';

folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_orig\';
outfolder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_computeBCs\';

folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_interpBCs_static_tail\';
outfolder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_interpBCs_static_tail\';

% folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_interpBCs_flipped_tail\';
% outfolder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_interpBCs_flipped_tail\';

folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_flipped_redo\';
outfolder = folder;


folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_orig_lateral_tail_angle\';
outfolder = folder;



static_tail = true;
 time_tail = 0;  %must match tail mesh file we are using for all time
flipped_tail = false;

static_tail = false;
flipped_tail = false;  % vertical tail angle
lateral_tail = true;  % horizontal tail angle

if ~exist(outfolder,'dir')
    mkdir(outfolder);
end


max_transverse_tries = 100;
% names = {'Body','Transverse','Tail'};
names = {'Body', 'Transverse', 'Tail', 'Wingtip'};
names_L = {'body', 'transverse', 'tail', 'wingtip'};

% names = { 'Transverse', 'Wingtip'};
% names_L = { 'transverse', 'wingtip'};

do_figure = false;
do_vid = false;

if do_vid
    
    profile = 'MPEG-4';
    profile = 'Motion JPEG AVI';
    %     vidh = VideoWriter('C:\Users\rudi\Desktop\RD\transverse_wingtip_vid',profile);
    vidh = VideoWriter('E:\Hull\dinoflagellate\transverse_wingtip_wide_vid',profile);
    
    vidh.FrameRate = 128/10; %50;
    vidh.Quality = 100; %1-100
    open(vidh);
    
end


store_constants = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear geom

geom.body.center = [10.7 * 0.825 / 2 + 10.7 / 2    0      0]';

geom.phase_speed = 2*pi*46;  %if lowest common period T = 4*pi s, covered in 2*pi rad, is 1/2 rad/sec for the entire beat cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geom.tail.radius = 0.15   * 1;
% geom.tail.radius = 0.9;


geom.tail.lambda = 22.2;
geom.tail.nlambda = 1.8531;
geom.tail.amp = [1.752*0.95  11.4 / 2];
% geom.tail.kE = [3 * 0.24      2 * 0.05 * 1.5 ];
geom.tail.kE = [3 * 0.24  * 0.5     2 * 0.05 * 1.5 ];


geom.tail.omega = 2 * pi * 46 ;   %rad/sec
% period = 2*pi rad / (2*pi*46  rad/sec) = 1/46 sec
% phase speed = 2*pi / period = 2*pi*46  rad/sec

% geom.tail.t_transition = 3.050645083561573;
% geom.tail.t_transition = 2.999558294702207;
geom.tail.t_transition = 3.000076881165709;

% geom.tail.translation = [8.326846940519667 0.000000000000004 14.100000000000000]';
% geom.tail.translation = [8.146345068104413 0.000000000000004 14.650000000000000]';
geom.tail.translation = [8.148177358404144 0.000000000000004 14.650000000000000]';

geom.tail.t_min = 0;
geom.tail.t_max = 2*pi * geom.tail.nlambda;  %t value where hemispherical end cap is centered
safety_factor = 0.1;  % t value at start of startign sph or end of ending sphere must be within (t_max - t_min)*safety_factor of t_min or t_max respectively

%geom.tail_angle = -25 *0  ;

geom.tail.rotation_angle = [-pi/2 ]';  %if multiple values, several consecutive rotations are done
% -90 deg vs + 90 deg in Salome because in Salome order of rotation is
% x then z, while here the z rotation is hardcoded into tail eqs and
% occurs first
geom.tail.rotation_pt  = [0 0 0]';
geom.tail.rotation_vec = [1 0 0]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geom.transverse.R = 13.45;  %orig before wingtip
geom.transverse.R = 13.45 + 0.5;  %shallow groove for wingtip

%geom.R = 15.7;  big groove I think
geom.transverse.d = 0;
%geom.w = 0.1;
geom.transverse.w = 0.6;  %4
geom.transverse.rf = 3.2 / 2;
geom.transverse.c = 18;
geom.transverse.omega = 2 * pi * 46;

geom.transverse.b = 3.7 ;  %wavelength of underlying helix - 0 for no helix, just a circle
geom.transverse.N_revs = 1;
%     geom.u_max = 2*pi*geom.N_revs;
% geom.transverse.u_max = 6.130414606345004; %no longer based on N_revs since sheet stops where groove ends, which Salome determines
geom.transverse.u_max = 6.126525207868735;

geom.transverse.v_max = 1;
%     geom.shift = 8.7;
geom.transverse.shift = 10.7 * 0.825 * 1.25;

%%%%%%%%%%%%%%%%%%%%%%
geom.transverse.wingtip_type = 'parallel_2';
geom.transverse.h_min = 0;  % always 0 for parallel wingtip
geom.transverse.h_max = 3;

% geom.transverse.h_min = -2;  %appears to be upper edge for normal wingtip
% geom.transverse.h_max = 2;  %appears to be lower edge for normal wingtip
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Inds Times Files SizeOK

for n = 1:length(names)
    name = names{n};
    
    files = dir([folder,name,'*']);
    sizes = [files.bytes];
    
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
    
    temp = sizes(sizes > 10E3);
    size_cutoffs = median(temp)*[ 1 - 0.5   ,  1 + 0.5 ];
    sizeOK = sizes >= size_cutoffs(1) & sizes <= size_cutoffs(2);
    SizeOK(:,n) = sizeOK(inds);
end

% if ~isequal(Times{1},Times{2},Times{3})
%     disp('problemo')
%     pause
% else
Time = Times{1};
% end





%%  get BCs for tail and then transverse for one phase angle at a time

for f = 1:length(Files{1})
    % fucked = 11 17 34 35 37 38 39 40 54 55 60 61 82 83 98 99 100
    
    if ~all(SizeOK(f,:))  %all files aren't present and good size
        disp(['Skipping step ',num2str(f),', file missing or size no good']);
        continue
    end
    
    disp(['Step ',num2str(f),' out of ',num2str(length(Files{1}))]);
    
    
      ind1 = strfind(Files{1}{f},'step_');
    if exist([outfolder,'Metadata_',Files{1}{f}(ind1:end-4),'.mat'],'file')
        disp('already done, skipping');
        continue
    end
    
    time = Time(f);
    clear Mesh Metadata_temp
    
    for i = 1:length(Files)
        [Mesh(i),Metadata_temp(i)] = load_mesh([folder,Files{i}{f}],[],[],'mesh');  %save rand vertex permutation for use later, so all future calculations use same exact mesh!
    end
    
    for i = 1:length(Files)
        Mesh(i).name = names_L{i};
        
        Mesh(i).orientation = [1 0 0; 0 1 0; 0 0 1;]';
        Mesh(i).refpoints = [0 0 0]';
    end
    

Mesh_orig = Mesh;

    
    
    ind = find(strcmp('tail',names_L));  %Mesh ind going with tail, if we loaded tail
    if ~isempty(ind)
        other_inds = setdiff(1:length(names_L), ind);  % inds for all other submeshes
        %tail indices might overlap body/transverse/wingtip since they're not part of same
        %mesh
        other_orig_elem = [];  other_orig_vert = [];
        for other_ind = other_inds
            other_orig_elem = [other_orig_elem; Mesh(other_ind).indices.orig.elem];
            other_orig_vert = [other_orig_vert; Mesh(other_ind).indices.orig.vert];
        end
        
        
        Mesh(ind).indices.orig.elem = Mesh(ind).indices.orig.elem + max(other_orig_elem);
        Mesh(ind).indices.orig.vert = Mesh(ind).indices.orig.vert + max(other_orig_vert);
        Mesh(ind).elems             = Mesh(ind).elems             + max(other_orig_vert);  %still using orig vert labels here, until renumber_Mesh below
    end
    

    
    
    [Mesh] = global_inds(Mesh);
    [Mesh] = renumber_Mesh(Mesh);
    
    if store_constants
        
        temp_input.performance.nthreads = input.performance.nthreads;
        
        temp_input.accuracy.integration_tol.area.abstol = 1E-9 ;
        temp_input.accuracy.integration_tol.area.reltol = 0.1;
        temp_input.accuracy.integration_tol.area.maxevals = Inf;
        
        temp_input.accuracy.integration_tol.centroid.abstol = 1E-9;
        temp_input.accuracy.integration_tol.centroid.reltol = 0.1;
        temp_input.accuracy.integration_tol.centroid.maxevals = Inf;
        
        temp_input.accuracy.integration_tol.volume.abstol = 1E-9;
        temp_input.accuracy.integration_tol.volume.reltol = 0.1;
        temp_input.accuracy.integration_tol.volume.maxevals = Inf;
        
        [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
    end
    
    if ismember('tail',names_L)
        [is_intersected] = submesh_intersections(Mesh, geom.tail.radius / 8 ,[], true, input.performance.nthreads);
        if is_intersected
            stopafra
        end
    end
    
    
    %% tail BCs
    if ismember('tail',names_L)
        tic;
        temp = {Mesh.name};
        ni = find(strcmp(temp,'tail'));
        clear vert_exact_tail
        t_sol = NaN(1,Mesh(ni).n_vert);
        dist = t_sol; theta = dist;  tail_vel = NaN(Mesh(ni).n_vert,3);  n_t = tail_vel;
        pts = tail_vel;  V_surfaces = pts;  V_us = pts;
        tail_errors = NaN(length(t_sol),2);
        disp('starting tail');
        parameters_tail = NaN(length(t_sol),2);
        parfor(vert_i = 1:length(t_sol), input.performance.nthreads)
%             disp(['tail vert ',num2str(vert_i)]);
            %  vert_i / Mesh(ni).n_vert
            vert = Mesh(ni).verts(vert_i,:);
            %   vert = Surface_vis(ni).verts(vert_i,:);
            if static_tail
               
                 distance = @(t) sqrt(sum((vert - tail_parameterized(t, time_tail, geom.tail)).^2)); 
            else
            distance = @(t) sqrt(sum((vert - tail_parameterized(t, time, geom.tail)).^2));  %distance between vert and tail centerline as function of t
            end
            t_sol(vert_i) = fminbnd(distance,geom.tail.t_min - (geom.tail.t_max - geom.tail.t_min)*(safety_factor)   ,   geom.tail.t_max + (geom.tail.t_max - geom.tail.t_min)*(safety_factor),optimset('maxfunevals',1E8,'maxiter',1E8,'tolx',1E-12)); %limit upper bound by t_max so points on hemisphere will fall at this parameter value
            
            
            dist(vert_i) = distance(t_sol(vert_i));
            
            %for spherical end cap, add vel of centerline pt to time der of vector
            %formed by rotating V_u in 3D around a new vector normal to both V_u and
            %V_surface (gotten with a cross product of the two) by an angle from the
            %dot product of the two
            if static_tail
                [pt, vel, der, V_u, radius] = tail_parameterized(t_sol(vert_i), time_tail, geom.tail);
            else
                [pt, vel, der, V_u, radius] = tail_parameterized(t_sol(vert_i), time, geom.tail);  % closest centerline pt and parameters of that centerline pt
            end
            pts(vert_i,:) = pt;
            % V_u is a more or less upward pointing vector that is normal
            % to the centerline (i.e. normal to der)
            V_u = V_u / sqrt(sum(V_u.^2));
            
            V_surface = vert - pt;  %vector from centerline pt to surface vert - depends on tail radius
            V_surface = V_surface / sqrt(sum(V_surface.^2));  %normalized direction only, doesn't depend on tail radius
            
            vert_exact_tail(vert_i,:) = real( pt + radius * V_surface );  % recalculate exact vert coord from closest pt on centerline and unit vector that should be normal to centerline pointed toward original vert
            % sometimes vel has a tiny imag component so take real part of
            % this too just to be safe
            
            
            V_surfaces(vert_i,:) = V_surface;
            V_us(vert_i,:) = V_u;
            
            
            
            %             theta(vert_i) = acos(dot(V_surface,V_u)) * (- sign(V_surface(3)));
            side = sign(V_surface(2));  %are we on the +y or -y side of the tail surface?
            if side == 0
                side = 1;  %if somehow we have a vert exact at the top or bottom, simply using 1 should work
            end
            theta(vert_i) = acos(dot(V_surface,V_u)) *  side;  %if on +y side, angle will be positive, and rotating V_u around der by theta will give V_surface direction
            % to obtain vert, take V_u from t_sol and rotate it around der
            % by theta and finally scale distance by local tail radius
            parameters_tail(vert_i,:) = [t_sol(vert_i) theta(vert_i)];
            
            %V_rot = cross(V_u, V_surface);
            if static_tail
                tail_vel(vert_i,:) = [0 0 0];
            else
                [n_t(vert_i,:)] = surface_vector_vel(t_sol(vert_i), theta(vert_i), time, geom.tail);  %time derivative of direction vector from centerline point to surface point
                
                %             tail_vel(vert_i,:) = vel + dist(vert_i) * n_t(vert_i,:);
                
                tail_vel(vert_i,:) = real( vel + radius * n_t(vert_i,:) );  %sometimes there is an extremely small imaginary component
            end
            
            if flipped_tail
                tail_vel(vert_i,:) = rotate_arbitrary_axis(tail_vel(vert_i,:),[0 0 0]', [0 1 0]', pi/2);  %rotate around origin since this is a velocity vector
            end
            
               if lateral_tail  % using 25 degrees since measured was 24 degrees
                tail_vel(vert_i,:) = rotate_arbitrary_axis(tail_vel(vert_i,:),[0 0 0]', [0 0 1]', -25*pi/180);  %rotate around origin since this is a velocity vector
            end
            
            if sqrt(sum(vel.^2)) < 1E-10  %vel mag is basically zero
                %  n_t(vert_i,:)
                % pause
            end
            %V_surface
            % V_u
            % V(i,:) = [0 1 0; -1 0 0; 0 0 1] * der(i,:)';  %rotate tangent vector 90 deg around z-axis (which is normal to movement plane) to obtain "upward" vector
            
            % [pt, vel, der] = change_ref_frame(pt, vel, der, geom.tail)
        end
        factor  = 0.99;  %to make sure we don't look at points on flat end surface at t = tmax
        
        tail_error = abs([min(dist(t_sol < geom.tail.t_max*factor)) max(dist(t_sol < geom.tail.t_max*factor))] - geom.tail.radius);  %an estimate of error = difference between distances and actual radius
        
        tail_errors(f,:) = tail_error;
        
        %         if max(tail_error) > 0.5
        %             disp('stopafra')
        %             pause
        %         end
        timings.tail = toc;
        
           if flipped_tail
                vert_exact_tail = rotate_arbitrary_axis(vert_exact_tail, geom.tail.translation, [0 1 0]', pi/2);
           end
           
              if lateral_tail
                vert_exact_tail = rotate_arbitrary_axis(vert_exact_tail, geom.tail.translation, [0 0 1]', -25*pi/180);
           end
            
        Mesh(ni).verts = vert_exact_tail;
        
        % make sure shared verts that were altered match exactly across all submeshes
        for nii = setdiff(1:length(Mesh),ni)
            [is_shared, inds] = ismember( Mesh(nii).indices.orig.vert , Mesh(ni).indices.orig.vert );
            
            Mesh(nii).verts( is_shared  , :) = Mesh(ni).verts( inds(is_shared), :);
        end
        
    end
    
    
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
    if ismember('transverse',names_L)
        tic;
        temp = {Mesh.name};
        ni = find(strcmp(temp,'transverse'));
        clear dist u_sol v_sol
        clear vert_exact_transverse
        u_sol = NaN(1,Mesh(ni).n_vert);
        %t_sol = NaN(1,Surface_vis(ni).n_vert);
        disp_probs = false;
        parameters_transverse = NaN(length(u_sol),2);
        dist = u_sol; v_sol = dist; transverse_vel = NaN(Mesh(ni).n_vert,3);
        dist_tol = 1E-1;
        tries = zeros(size(dist));
        parameter_tol = 1E-2;
        %     tol = 0.1;
        param_err = NaN(Mesh(ni).n_vert,2);
        transverse_dist_error = NaN(length(u_sol),1);
        transverse_param_error = NaN(length(u_sol),2);
        disp('starting transverse');
        meshbad = false;
        parfor (vert_i = 1:length(u_sol), input.performance.nthreads)
            %                         for vert_i = 1:length(u_sol)
            
            
            
           % disp(['transverse vert ',num2str(vert_i)]);
            %  vert_i / Mesh(ni).n_vert
            vert = Mesh(ni).verts(vert_i,:);
            %   vert = Surface_vis(ni).verts(vert_i,:);
            
            distance = @(x) sqrt(sum((vert' - transverse_parameterized(x(1),x(2), time, geom.transverse)).^2));
            % answer = fmincon(distance,[0,0],[],[],[],[],[0 0],[geom.transverse.u_max geom.transverse.v_max],[],optimset('maxfunevals',1E20,'maxiter',1E20,'tolx',1E-13,'tolfun',1E-13)); %limit upper bound by t_max so points on hemisphere will fall at this parameter value
            dist(vert_i) = Inf; u_sol(vert_i) = Inf;  v_sol(vert_i) = Inf;
            
            trycount = 0;
            while dist(vert_i) > dist_tol || u_sol(vert_i) < 0 - parameter_tol || v_sol(vert_i) < 0 - parameter_tol || u_sol(vert_i) > geom.transverse.u_max + parameter_tol || v_sol(vert_i) > geom.transverse.v_max + parameter_tol
                %                 answer = fminsearch(distance,[geom.transverse.u_max*rand, geom.transverse.v_max*rand],optimset('maxfunevals',1E2,'maxiter',1E3,'tolx',1E-6,'tolfun',1E-13))
                trycount = trycount + 1;
                answer = fmincon(distance, [geom.transverse.u_max*rand, geom.transverse.v_max*rand], [], [], [], [], [0 0]', [geom.transverse.u_max geom.transverse.v_max]', [],optimoptions('fmincon','display','off'));
                
                
                u_sol(vert_i) = answer(1);  v_sol(vert_i) = answer(2);
                parameters_transverse(vert_i,:) = [u_sol(vert_i) v_sol(vert_i)];
                
                
                vert_exact_transverse(vert_i,:) = transverse_parameterized(answer(1),answer(2), time, geom.transverse);
                
                
                dist(vert_i) = distance(answer);
%                 dist(vert_i)
                if u_sol(vert_i) < 0
                    u_err = 0 - u_sol(vert_i);
                elseif u_sol(vert_i) > geom.transverse.u_max
                    u_err = geom.transverse.u_max - u_sol(vert_i);
                else
                    u_err = 0;
                end
                if v_sol(vert_i) < 0
                    v_err = 0 - v_sol(vert_i);
                elseif v_sol(vert_i) > geom.transverse.v_max
                    v_err = geom.transverse.v_max - v_sol(vert_i);
                else
                    v_err = 0;
                end
                
                param_err(vert_i,:) = [u_err v_err];
                
                [~,transverse_vel(vert_i,:)] = transverse_parameterized(u_sol(vert_i),v_sol(vert_i), time, geom.transverse);
                tries(vert_i) =  tries(vert_i) + 1;
                %   [u_sol(vert_i)  v_sol(vert_i)]
                % dist(vert_i)
                if   dist(vert_i) > dist_tol && disp_probs
                    disp(['dist > dist_tol, dist = ',num2str(dist(vert_i))]);
                end
                
                if u_sol(vert_i) < 0 - parameter_tol  && disp_probs
                    disp(['u_sol < 0, u_sol = ',num2str(u_sol(vert_i))]);
                end
                
                if v_sol(vert_i) < 0 - parameter_tol  && disp_probs
                    disp(['v_sol < 0, v_sol = ',num2str(v_sol(vert_i))]);
                end
                
                if u_sol(vert_i) > geom.transverse.u_max + parameter_tol   && disp_probs
                    disp(['u_sol > u_max, u_sol - u_max = ',num2str(u_sol(vert_i) - geom.transverse.u_max)]);
                end
                
                if v_sol(vert_i) > geom.transverse.v_max + parameter_tol  && disp_probs
                    disp(['v_sol > v_max, v_sol - v_max = ',num2str(v_sol(vert_i) - geom.transverse.v_max)]);
                end
                
                if trycount > max_transverse_tries
                    meshbad = meshbad | true;
                    break;
                end
                
            end  %while
            
            
            
            
        end  %parfor verts
        
        
        if meshbad
            continue;
        end
        
        transverse_dist_error(f) = max(dist);
        transverse_param_error(f,:) = max(param_err);
        
        %         if transverse_dist_error(f) > 0.1
        %             disp('stopafra')
        %             pause
        %         end
        timings.transverse = toc;
        
        Mesh(ni).verts = vert_exact_transverse;
        
        for nii = setdiff(1:length(Mesh),ni)
            [is_shared, inds] = ismember( Mesh(nii).indices.orig.vert , Mesh(ni).indices.orig.vert );
            
            Mesh(nii).verts( is_shared  , :) = Mesh(ni).verts( inds(is_shared), :);
        end
    end
    
    
    %%
    if ismember('wingtip',names_L)
        tic;
        temp = {Mesh.name};
        ni = find(strcmp(temp,'wingtip'));
        clear dist u_sol h_sol
        u_sol = NaN(1,Mesh(ni).n_vert);
        %t_sol = NaN(1,Surface_vis(ni).n_vert);
        disp_probs = false;
        
        dist = u_sol; h_sol = dist; wingtip_vel = NaN(Mesh(ni).n_vert,3);
        parameters_wingtip = NaN(length(u_sol),2);
        %dist_tol = 5E-2;
        %         dist_tol = 1E-1;  %sometimes get 0.4 off...
        
        dist_tols = [0.05 0.1 0.15 0.2];  %parallel_2
          dist_tols = [0.01  0.05 0.1  0.25 ];  %parallel_2
        %          dist_tol = 0.04;  %parallel_1
        %          dist_tol = 0.2;
        max_tries = 20;
        tries = zeros(size(dist));
        parameter_tol = 5E-2;
        %     tol = 0.1;
        clear vert_exact_wingtip
        param_err = NaN(Mesh(ni).n_vert,2);
        wingtip_dist_error = NaN(length(u_sol),1);
        wingtip_param_error = NaN(length(u_sol),2);
        disp('starting wingtip');
        meshbad = false;
        parfor (vert_i = 1:length(u_sol), input.performance.nthreads)
            answer = [];
            %             if dist(vert_i) < 0.03
            %                 continue
            %             end
            
                   %   for vert_i = 1:length(u_sol)
%             disp(['wingtip vert ',num2str(vert_i)]);
            %  vert_i / Mesh(ni).n_vert
            vert = Mesh(ni).verts(vert_i,:);
            %   vert = Surface_vis(ni).verts(vert_i,:);
            
            distance = @(x) sqrt(sum((vert' - transverse_hairs_parameterized_mexed(x(1),x(2), time, geom.transverse)).^2));
            % answer = fmincon(distance,[0,0],[],[],[],[],[0 0],[geom.transverse.u_max geom.transverse.v_max],[],optimset('maxfunevals',1E20,'maxiter',1E20,'tolx',1E-13,'tolfun',1E-13)); %limit upper bound by t_max so points on hemisphere will fall at this parameter value
            dist(vert_i) = Inf; u_sol(vert_i) = Inf;  h_sol(vert_i) = Inf;
            
            still_bad = true;
            
            for d = 1:length(dist_tols)
                if ~still_bad  %success
                    break
                end
                
                dist_tol = dist_tols(d);
                
                
                
                tc = 0;  still_bad = false;
                while dist(vert_i) > dist_tol || u_sol(vert_i) < 0 - parameter_tol || h_sol(vert_i) < geom.transverse.h_min - parameter_tol || u_sol(vert_i) > geom.transverse.u_max + parameter_tol || h_sol(vert_i) > geom.transverse.h_max + parameter_tol
                    % answer = fminsearch(distance,[geom.transverse.u_max*rand, geom.transverse.h_min + (geom.transverse.h_max - geom.transverse.h_min)*rand],optimset('maxfunevals',1E20,'maxiter',1E4,'tolx',1E-13,'tolfun',1E-13));
                    tc = tc + 1;
                    
                    if tc > max_tries
                        still_bad = true;
                        
                        break
                    end
                    
                    %                 answer = fmincon(distance, [geom.transverse.u_max*rand, geom.transverse.h_min + (geom.transverse.h_max - geom.transverse.h_min)*rand], [],[],[],[],[0 geom.transverse.h_min]',[geom.transverse.u_max geom.transverse.h_max]',[],optimoptions('fmincon','display','off'))
                    try
                        answer = patternsearch(distance, [geom.transverse.u_max*rand, geom.transverse.h_min + (geom.transverse.h_max - geom.transverse.h_min)*rand], [],[],[],[],[0 geom.transverse.h_min]',[geom.transverse.u_max geom.transverse.h_max]',[],optimoptions('patternsearch','display','off','PollOrderAlgorithm', 'random', 'UseCompletePoll', false,'AccelerateMesh',true,'cache','on'));
%                         answer = fmincon(      distance, [geom.transverse.u_max*rand, geom.transverse.h_min + (geom.transverse.h_max - geom.transverse.h_min)*rand], [],[],[],[],[0 geom.transverse.h_min]',[geom.transverse.u_max geom.transverse.h_max]',[],optimoptions('fmincon','display','off'));
                 
                    catch
                        answer = NaN(1,2);
                    end
                    
                    
                    temp =  transverse_hairs_parameterized_mexed(answer(1),answer(2), time, geom.transverse);
                   
                     vert_exact_wingtip(vert_i,:) = temp;
                     
                    if any(isnan(  temp ) );
                        stoparoo
                    end
                    %                 [ vert_i   vert_i / Mesh(ni).n_vert]
                    u_sol(vert_i) = answer(1);  h_sol(vert_i) = answer(2);
                    
                    parameters_wingtip(vert_i,:) = [u_sol(vert_i) h_sol(vert_i)];
                    
                    dist(vert_i) = distance(answer);
                    
                   % [answer dist(vert_i)]
                    
                    if u_sol(vert_i) < 0
                        u_err = 0 - u_sol(vert_i);
                    elseif u_sol(vert_i) > geom.transverse.u_max
                        u_err = geom.transverse.u_max - u_sol(vert_i);
                    else
                        u_err = 0;
                    end
                    
                    if h_sol(vert_i) < geom.transverse.h_min
                        h_err = geom.transverse.h_min - h_sol(vert_i);
                    elseif h_sol(vert_i) > geom.transverse.h_max
                        h_err = geom.transverse.h_max - h_sol(vert_i);
                    else
                        h_err = 0;
                    end
                    
                    param_err(vert_i,:) = [u_err h_err];
                    
                    [~,wingtip_vel(vert_i,:)] = transverse_hairs_parameterized_mexed(u_sol(vert_i),h_sol(vert_i), time, geom.transverse);
                    tries(vert_i) =  tries(vert_i) + 1;
                    %   [u_sol(vert_i)  v_sol(vert_i)]
                    % dist(vert_i)
                    if   dist(vert_i) > dist_tol && disp_probs
                        disp(['dist > dist_tol, dist = ',num2str(dist(vert_i))]);
                    end
                    
                    if u_sol(vert_i) < 0 - parameter_tol  && disp_probs
                        disp(['u_sol < 0, u_sol = ',num2str(u_sol(vert_i))]);
                    end
                    
                    if h_sol(vert_i) < geom.transverse.h_min - parameter_tol  && disp_probs
                        disp(['h_sol < ',num2str(geom.transverse.h_min),', h_sol = ',num2str(h_sol(vert_i))]);
                    end
                    
                    if u_sol(vert_i) > geom.transverse.u_max + parameter_tol   && disp_probs
                        disp(['u_sol > u_max, u_sol - u_max = ',num2str(u_sol(vert_i) - geom.transverse.u_max)]);
                    end
                    
                    if h_sol(vert_i) > geom.transverse.h_max + parameter_tol  && disp_probs
                        disp(['h_sol > ',num2str(geom.transverse.h_max),', h_sol = ',num2str(h_sol(vert_i))]);
                    end
                    
                    
                    
                end  %while
                
            end  %for dist_tols
            
            if d == length(dist_tols) && still_bad
                meshbad = meshbad | true;
            end
            
            
        end  %parfor verts
        
        if meshbad
           
            continue;
        end
        %%
        wingtip_dist_error(f) = max(dist);
        wingtip_param_error(f,:) = max(param_err);
        
        %         if wingtip_dist_error(f) > dist_tol
        %             disp('stopafra')
        %             pause
        %         end
        timings.wingtip = toc;
        
        Mesh(ni).verts = vert_exact_wingtip;
        
        for nii = setdiff(1:length(Mesh),ni)
            [is_shared, inds] = ismember( Mesh(nii).indices.orig.vert , Mesh(ni).indices.orig.vert );
            
            Mesh(nii).verts( is_shared  , :) = Mesh(ni).verts( inds(is_shared), :);
        end
        
    end  %wingtip
    
    
if any(isnan(Mesh(4).verts(:)))
    stoperdoodle
end
  
%         stoperdoodle
    %%
    clear Metadata
    Metadata.geom = geom;
    Metadata.time = time;
    
    for i = 1:length(names_L)  %ensures order within Metadata matches order in which submeshes were loaded (which is order of names_L)
        name = names_L{i};
        switch name
            case 'body'
                Metadata.body.BCs = zeros(Mesh(i).n_vert,3);
                
            case 'transverse'
                Metadata.transverse.BCs = transverse_vel;
                Metadata.transverse.parameters = parameters_transverse;
            case 'tail'
                Metadata.tail.BCs = tail_vel;
                Metadata.tail.parameters = parameters_tail;
            case 'wingtip'
                Metadata.wingtip.BCs = wingtip_vel;
                Metadata.wingtip.parameters = parameters_wingtip;
        end
        Metadata.(name).indices.orig.vert = Mesh(i).indices.orig.vert;
        %         Metadata.(name).rand_inds = Metadata_temp(i).mesh.rand_inds;
      %  Metadata.(name).rand_inds = [1:Mesh(i).n_vert];
    end
    
    ind1 = strfind(Files{1}{f},'step_');
    save([outfolder,'Metadata_',Files{1}{f}(ind1:end-4),'.mat'],'Metadata');
    
    
    for i = 1:length(names_L)
        
        
        filename = [outfolder, names{i},'_',Files{1}{f}(ind1:end)];
        if exist(filename,'file')
            delete(filename);
        end
        dlmwrite(filename,[Mesh(i).n_vert Mesh(i).n_elem],'-append','delimiter',' ','precision',16)
        dlmwrite(filename,[Mesh_orig(i).indices.orig.vert    Mesh(i).verts],'-append','delimiter',' ','precision',16)
        dlmwrite(filename,[Mesh_orig(i).indices.orig.elem      repmat(206,Mesh(i).n_elem,1)     Mesh_orig(i).elems],'-append','delimiter',' ','precision',16)
        % we shifted up orig vert and elem inds as well as actual elems and
        % need to use those here instead of current Mesh.elems since this
        % has been renumbered to start at 1 for each submesh
    end
    
    
    
    
    
    
    
    %%
    if do_figure
        figure(488)
        clf
        %     if ismember('body',names_L) && ismember('tail',names_L) && ismember('transverse',names_L)
        %         refines = [3 2 1];
        %     elseif ismember('transverse',names_L) && ismember('wingtip', names_L)
        %         refines = [2 2];
        %     end
        [s,e] = plot_mesh(Mesh, [2 2 2 2]);
        hold on
        
        %%
        
        %     try, delete(q1), end;   try, delete(q2), end;  try, delete(q3), end;
        %
        %     if ismember('tail',names_L)
        %         inc = 1;  ni = find(strcmp('tail',names_L));
        %         q1 = quiver3(Mesh(ni).verts(1:inc:end,1),Mesh(ni).verts(1:inc:end,2),Mesh(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3), 1);
        %         %q = quiver3(Surface_vis(ni).verts(1:inc:end,1),Surface_vis(ni).verts(1:inc:end,2),Surface_vis(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3));
        %         set(q1,'color','k');
        %     end
        %     if ismember('transverse',names_L)
        %         inc = 1;  ni = find(strcmp('transverse',names_L));
        %         q2 = quiver3(Mesh(ni).verts(1:inc:end,1),Mesh(ni).verts(1:inc:end,2),Mesh(ni).verts(1:inc:end,3),transverse_vel(1:inc:end,1),transverse_vel(1:inc:end,2),transverse_vel(1:inc:end,3), 1);
        %         %q = quiver3(Surface_vis(ni).verts(1:inc:end,1),Surface_vis(ni).verts(1:inc:end,2),Surface_vis(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3));
        %         set(q2,'color','k');
        %     end
        %     if ismember('wingtip',names_L)
        %         inc = 1;  ni = find(strcmp('wingtip',names_L));
        %         q3 = quiver3(Mesh(ni).verts(1:inc:end,1),Mesh(ni).verts(1:inc:end,2),Mesh(ni).verts(1:inc:end,3),wingtip_vel(1:inc:end,1),wingtip_vel(1:inc:end,2),wingtip_vel(1:inc:end,3), 1);
        %         %q = quiver3(Surface_vis(ni).verts(1:inc:end,1),Surface_vis(ni).verts(1:inc:end,2),Surface_vis(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3));
        %         set(q3,'color','k');
        %     end
        
        
        %       inc = 1;
        %         q4 = quiver3([Mesh(1).verts(1:inc:end,1); Mesh(2).verts(1:inc:end,1)],[Mesh(1).verts(1:inc:end,2); Mesh(2).verts(1:inc:end,2)],[Mesh(1).verts(1:inc:end,3); Mesh(2).verts(1:inc:end,3)],[transverse_vel(1:inc:end,1); wingtip_vel(1:inc:end,1)],[transverse_vel(1:inc:end,2); wingtip_vel(1:inc:end,2)],[transverse_vel(1:inc:end,3); wingtip_vel(1:inc:end,3)], 1);
        %         %q = quiver3(Surface_vis(ni).verts(1:inc:end,1),Surface_vis(ni).verts(1:inc:end,2),Surface_vis(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3));
        %         set(q4,'color','k');
        
        transverse_ind = find(strcmp({Mesh.name},'transverse'));
        tail_ind = find(strcmp({Mesh.name},'tail'));
       wingtip_ind = find(strcmp({Mesh.name},'wingtip'));
        try, delete(q4); end
        inc = 1;  factor = 2;
        q4 = quiver3([ Mesh(transverse_ind).verts(1:inc:end,1); Mesh(tail_ind).verts(1:inc:end,1); Mesh(wingtip_ind).verts(1:inc:end,1)],[ Mesh(transverse_ind).verts(1:inc:end,2); Mesh(tail_ind).verts(1:inc:end,2); Mesh(wingtip_ind).verts(1:inc:end,2)],[Mesh(transverse_ind).verts(1:inc:end,3); Mesh(tail_ind).verts(1:inc:end,3); Mesh(wingtip_ind).verts(1:inc:end,3)],[ transverse_vel(1:inc:end,1); tail_vel(1:inc:end,1); wingtip_vel(1:inc:end,1)],[ transverse_vel(1:inc:end,2); tail_vel(1:inc:end,2); wingtip_vel(1:inc:end,2)],[transverse_vel(1:inc:end,3); tail_vel(1:inc:end,3); wingtip_vel(1:inc:end,3)], factor);
        %q = quiver3(Surface_vis(ni).verts(1:inc:end,1),Surface_vis(ni).verts(1:inc:end,2),Surface_vis(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3));
        set(q4,'color','k');
        
        %     xlim([-36       26.717]);  ylim([   -23.663       23.664]);  zlim([  -20.705       23.301]);
        xlim([0       20]);  ylim([   -23.663       23.664]);  zlim([  -20.705       23.301]);
        
        %     set(gca,'view',[  -24.321       47.549]);
        set(gca,'view',[  -80      78.8]);
        
        
        set(gcf,'position',[  582          31        1292         974]);
        set(gca,'cameraViewAngle',4.0468);
        set(gca,'cameraPosition',[-128.15      -267.48       327.41]);
        set(gca,'cameraTarget',[-5.8668       3.0753       2.8308]);
        
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
end  %time steps e.g. filenames


if do_vid
    close(vidh)
end