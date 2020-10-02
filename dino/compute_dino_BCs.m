
min_free_mem = 10; % min free memory (GB) before deleting and restarting parpool to fix memory leak
input.performance.nthreads = 20;


progress_monitors = false;


time_tail = 0;  %must match tail mesh file we are using for all time


static_tail = false;
flipped_tail = false;  % vertical tail angle
lateral_tail = false;  % horizontal tail angle
centered_tail = false;
% tail_as_sheet = true;
tail_as_sheet = false;
fat_tail = false;  % radius 0.15 vs 0.75 micron

tail_shift = 0;  % how many microns to shift tail by (in distal direction, away from transverse groove and sheets)
% needed for fat cylindrical tail to not hit hairs at around step 64 / 128;
%    1 seems to work, however we don't even need this for centered tail
%    case anyway

actively_flipping_tail = false;

% names = {'Body','Transverse','Tail'};
names = {'Body', 'Transverse', 'Tail', 'Coplanar_Hairs' , 'Normal_Top_Hairs' , 'Normal_Bottom_Hairs'};
% names_L = {'body', 'transverse', 'tail', 'wingtip' , 'coplanar-hairs' , 'normal-top-hairs' , 'normal-bottom-hairs'};

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
% dino_geom_parameters_huge_groove_uber_hairs_longtail;
% dino_geom_parameters_huge_groove_uber_hairs_3p5_1p5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if actively_flipping_tail
    %     load('C:\Users\rudi\Desktop\RD\dino turning\tail_flipping_spline.mat'); % tail flips out from body over a beat cycle (loads tail_spline variable)
    load('C:\Users\rudi\Desktop\RD\dino turning\tail_unflipping_spline.mat'); % tail unflips back to regular over a beat cycle
    temp = 0; % keep final step e.g. 128 in since it isn't the same as first step
else
    temp = 1; % remove final step since it's the same as the first
end

step_list = 0:(n_phase_pts - temp);

total_period = 1 / 46.0 ;
time_list = linspace(0, total_period, n_phase_pts + 1)  ;
time_list = time_list(1:end - temp);  % leave off last value which is same as first
phase_list = linspace(0, 2*pi, n_phase_pts + 1);   % *always* 2 pi rad in a cycle
phase_list = phase_list(1:end - temp);  % leave off last value which is same as first





clear Steps Phases Times Files SizeOK

for n = 1:length(names)
    name = names{n};
    
    files = dir([folder_original,name,'*']);
    sizes = [files.bytes];
    
    files = {files.name};
    %files = files(3:end);
    if isempty(files)
        return;
    end
    
    clear times
    for i = 1:length(files)
        
        ind = strfind(files{i},'step_');
        temp = sscanf(files{i}(ind+5:end),'%f_time_%f_phase_%f.dat'); %[step, time, phase]
        times(i) = temp(2);
        indices(i) = temp(1);
        phases(i) =temp(3);
    end
    
    [~,inds] = sort(times);
    %     Inds{n} = inds;
    Files{n} = files(inds);
    Times{n} = times(inds);
    Steps{n} = indices(inds);
    Phases{n} = phases(inds);
    
    %     if ~tail_mesh_only
    temp = sizes(sizes > 10E3);
    size_cutoffs = median(temp)*[ 1 - 0.5   ,  1 + 0.5 ];
    sizeOK = sizes >= size_cutoffs(1) & sizes <= size_cutoffs(2);
    %     SizeOK(:,n) = sizeOK(inds);
    SizeOK{n} = sizeOK(inds);
    %     end
end

% if ~isequal(Times{1},Times{2},Times{3})
%     disp('problemo')
%     pause
% else
Time = Times{1};
% end



last_solution.Coplanar_Hairs = [];  last_solution.Normal_Top_Hairs = [];  last_solution.Normal_Bottom_Hairs = [];

%%  get BCs for tail and then transverse for one phase angle at a time
% if tail_mesh_only
%     file_ind = 3;
% else
file_ind = 1; % Body meshes control which indices get done
% end
f_list = 1:length(Files{file_ind});
%  f_list = [ f_list(randperm(length(f_list)))];
% f_list = f_list(1);
for f = f_list
    
    
    [~,temp] = system('wmic OS get FreePhysicalMemory /Value'); % memory status according to Windows
    free_mem = str2num( temp( regexp(temp,'\d')  ) ) / 1E6;  % free memory in GB
    if free_mem < min_free_mem
        delete(gcp);
        pause(5);
        parpool(input.performance.nthreads);
        pause(2);
    end
    
    %     if tail_mesh_only
    %     temp = sscanf(Files{file_ind}{f},'Tail_step_%f_time_%f_phase_%f.dat'); %[mesh step, time, phase]
    %     else
    temp = sscanf(Files{file_ind}{f},'Body_step_%f_time_%f_phase_%f.dat'); %[step, time, phase]
    %     end
    
    step = temp(1); phase = phase_list( (step == step_list) );
    clear local_index
    sizes_OK = NaN(length(names),1);
    for n = 1:length(names)
        [~, local_index(n)] = ismember(step,Steps{n});
        if local_index(n) ~= 0
            sizes_OK(n) = SizeOK{n}(local_index(n));
        end
    end
    
    
    
    if any(isnan(sizes_OK)) || ~all(sizes_OK)  %all files aren't present and good size
        %         disp(['Skipping step ',num2str(f),', file missing or size no good']);
        continue
    end
    
    %     disp(['Index ',num2str(f),' out of ',num2str(length(Files{1}))]);
    
    
    ind1 = strfind(Files{file_ind}{f},'step_');
    if exist([folder_processed,'Metadata_',Files{file_ind}{f}(ind1:end-4),'.mat'],'file')
        %         disp('already done, skipping');
        continue
    end
    
    parameters_file = [folder_original , 'Parameters_step_',num2str(temp(1)), '_time_',num2str(temp(2),'%.15f'),'_phase_',num2str(temp(3),'%.15f'),'.txt'];
    
    %     parameters_file = 'C:/Users/rudi/Desktop/RD/dino_mesh_automation/phases_hairs_4_1/original/Parameters_step_0_time_0.000000000000000_phase_0.000000000000000.txt';
    if ~exist(parameters_file,'file')
        disp('parameters_file doesn''t exist; trying to use file from any step');
        temp = dir([folder_original,'Parameters*.txt']);  temp = {temp.name};
        if ~isempty(temp)
            parameters_file = [folder_original, temp{1}];
        else
            error('no parameters files found');
        end
    end
    dino_geom_parameters_from_file;
    
    
    time = Time(f);
    clear Mesh Metadata_temp
    
    for i = 1:length(Files)
        [Mesh(i),Metadata_temp(i)] = load_mesh([folder_original,Files{i}{local_index(i)}],[],[],'mesh');  %save rand vertex permutation for use later, so all future calculations use same exact mesh!
    end
    
    for i = 1:length(Files)
        Mesh(i).name = names{i};
        
        Mesh(i).orientation = [1 0 0; 0 1 0; 0 0 1;]';
        Mesh(i).refpoints = [0 0 0]';
    end
    
    
    Mesh_orig = Mesh;
    
    
    
    ind = find(strcmp('Tail',names));  %Mesh ind going with tail, if we loaded tail
    if ~isempty(ind)
        other_inds = setdiff(1:length(names), ind);  % inds for all other submeshes
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
        clear temp_input
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
    
    
    meshbad = false;
    
    %% tail BCs
    if ismember('Tail',names)
        tic;
        temp = {Mesh.name};
        ni = find(strcmp(temp,'Tail'));
        clear vert_exact_tail
        t_sol = NaN(1,Mesh(ni).n_vert);
        dist = t_sol; theta = dist;  tail_vel = NaN(Mesh(ni).n_vert,3);  n_t = tail_vel;
        pts = tail_vel;  V_surfaces = pts;  V_us = pts;
        tail_errors = NaN(length(t_sol),2);
        disp('starting tail');
        parameters_tail = NaN(length(t_sol),2);
        parfor(vert_i = 1:length(t_sol), input.performance.nthreads)
            %  disp(['tail vert ',num2str(vert_i)]);
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
            
            if tail_as_sheet
                vert_exact_tail(vert_i,:) = real( [pt(1) vert(2) pt(3)]);  %trust across-sheet coordinate of vert, but a point on the sheet must have the same x and z coords at closest pt on centerline
            else
                vert_exact_tail(vert_i,:) = real( pt + radius * V_surface );  % recalculate exact vert coord from closest pt on centerline and unit vector that should be normal to centerline pointed toward original vert
                % sometimes vel has a tiny imag component so take real part of
                % this too just to be safe
            end
            
            V_surfaces(vert_i,:) = V_surface;
            V_us(vert_i,:) = V_u;
            
            
            if ~tail_as_sheet
                %             theta(vert_i) = acos(dot(V_surface,V_u)) * (- sign(V_surface(3)));
                side = sign(V_surface(2));  %are we on the +y or -y side of the tail surface?
                if side == 0
                    side = 1;  %if somehow we have a vert exact at the top or bottom, simply using 1 should work
                end
                theta(vert_i) = acos(dot(V_surface,V_u)) *  side;  %if on +y side, angle will be positive, and rotating V_u around der by theta will give V_surface direction
                % to obtain vert, take V_u from t_sol and rotate it around der
                % by theta and finally scale distance by local tail radius
                parameters_tail(vert_i,:) = [t_sol(vert_i) theta(vert_i)];
            else
                radius = vert(2) - pt(2);  %distance along sheet between centerline pt and vert, like a local hair length parameter
                parameters_tail(vert_i,:) = [t_sol(vert_i) radius];
            end
            
            %V_rot = cross(V_u, V_surface);
            if static_tail
                tail_vel(vert_i,:) = [0 0 0];
            else
                if ~tail_as_sheet
                    [n_t(vert_i,:)] = surface_vector_vel(t_sol(vert_i), theta(vert_i), time, geom.tail);  %time derivative of direction vector from centerline point to surface point
                    tail_vel(vert_i,:) = real( vel + radius * n_t(vert_i,:) );  %sometimes there is an extremely small imaginary component
                else
                    tail_vel(vert_i,:) = real( vel  );  %velocity is constant across sheet
                end
            end
            
            if flipped_tail
                tail_vel(vert_i,:) = rotate_arbitrary_axis(tail_vel(vert_i,:),[0 0 0]', [0 1 0]', pi/4);  %rotate around origin since this is a velocity vector
            end
            
            if actively_flipping_tail
                tail_vel(vert_i,:) = rotate_arbitrary_axis(tail_vel(vert_i,:),[0 0 0]', [0 1 0]', slmeval(phase, tail_spline));  %rotate around origin since this is a velocity vector
            end
            
            if lateral_tail  % using 25 degrees since measured was 24 degrees
                tail_vel(vert_i,:) = rotate_arbitrary_axis(tail_vel(vert_i,:),[0 0 0]', [0 0 1]', -25*pi/180);  %rotate around origin since this is a velocity vector
            end
            
            
            %             if sqrt(sum(vel.^2)) < 1E-10  %vel mag is basically zero
            %                 %  n_t(vert_i,:)
            %                 % pause
            %             end
            %             %V_surface
            % V_u
            % V(i,:) = [0 1 0; -1 0 0; 0 0 1] * der(i,:)';  %rotate tangent vector 90 deg around z-axis (which is normal to movement plane) to obtain "upward" vector
            
            % [pt, vel, der] = change_ref_frame(pt, vel, der, geom.tail)
        end
        factor  = 0.99;  %to make sure we don't look at points on flat end surface at t = tmax
        
        if ~tail_as_sheet
            tail_error = abs([min(dist(t_sol < geom.tail.t_max*factor)) max(dist(t_sol < geom.tail.t_max*factor))] - geom.tail.radius);  %an estimate of error = difference between distances and actual radius
            
            tail_errors(f,:) = tail_error;
        else
            tail_errors(f,:) = NaN;
        end
        %         if max(tail_error) > 0.5
        %             disp('stopafra')
        %             pause
        %         end
        timings.tail = toc;
        
        if flipped_tail
            vert_exact_tail = rotate_arbitrary_axis(vert_exact_tail, geom.tail.translation, [0 1 0]', pi/4);
        end
        
        if actively_flipping_tail
            vert_exact_tail = rotate_arbitrary_axis(vert_exact_tail, geom.tail.translation, [0 1 0]',  slmeval(phase, tail_spline));
            
        end
        
        if lateral_tail
            vert_exact_tail = rotate_arbitrary_axis(vert_exact_tail, geom.tail.translation, [0 0 1]', -25*pi/180);
        end
        
        if centered_tail
            vert_exact_tail  =  vert_exact_tail  + repmat(geom.tail.shift',size(vert_exact_tail,1),1);
        end
        
        vert_exact_tail = vert_exact_tail + repmat([-tail_shift 0 0],size(vert_exact_tail,1),1);
        
        
        Mesh(ni).verts = vert_exact_tail;
        
        % make sure shared verts that were altered match exactly across all submeshes
        for nii = setdiff(1:length(Mesh),ni)
            [is_shared, inds] = ismember( Mesh(nii).indices.orig.vert , Mesh(ni).indices.orig.vert );
            
            Mesh(nii).verts( is_shared  , :) = Mesh(ni).verts( inds(is_shared), :);
        end
        
        
        
        
        
        % need to test here, after we may have moved or rotated tail
        % outside Salome
        [is_intersected] = submesh_intersections(Mesh, geom.tail.radius / 8 ,[], true, input.performance.nthreads);
        if is_intersected
            error('Tail submesh intersects some other submesh');
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
    %
    %             u = linspace(0,geom.transverse.u_max,1000);
    %             v = linspace(0,geom.transverse.v_max,50);
    %
    %             v = linspace(geom.hairs.Coplanar_Hairs.h_min,geom.hairs.Coplanar_Hairs.h_max,100);
    %
    %
    %             c = 0;
    %             pts = NaN(length(u)*length(v),3);  U = NaN(length(u)*length(v),1); V = U;
    %             for i = 1:length(u)
    %                 i/length(u)
    %
    %                 for j = 1:length(v)
    %                     c = c+1;
    %                     U(c) = u(i);  V(c) = v(j);
    %
    % %                     pts(c,:) = transverse_parameterized(u(i),v(j), time, geom.transverse)';
    %                     [pts(c,:)] = transverse_hairs_parameterized(u(i),v(j), time, geom.transverse , 'Coplanar_Hairs');
    %                 end
    %             end
    %
    %             hold on
    %             ph = plot3(pts(:,1),pts(:,2),pts(:,3),'k.');
    
    %% show closest points on theoretical surface to each vertex of actual mesh
    %         u = u_sol;  v = v_sol;  %u and v had better be same size since they match with verts
    %
    %         %c = 0;
    %         pts_fit = NaN(length(u),3);  U = NaN(length(u),1); V = U;
    %         parfor i = 1:length(u)
    %             %  i/length(u)
    %
    %             % c = c+1;
    %             %             U(c) = u(i);  V(c) = v(j);
    %
    %             pts_fit(i,:) = transverse_parameterized(u(i),v(i), time, geom.transverse)';
    %
    %         end
    %
    %         hold on
    %         ph = plot3(pts_fit(:,1),pts_fit(:,2),pts_fit(:,3),'k.');
    
    %%
    %     clear th
    %     for i = 1:size(pts_fit,1)
    %         th(i) = text(pts_fit(i,1),pts_fit(i,2),pts_fit(i,3),num2str(v_sol(i)),'fontsize',12);
    %     end
    %%
    if ismember('Transverse',names)
        tic;
        temp = {Mesh.name};
        ni = find(strcmp(temp,'Transverse'));
        clear dist u_sol v_sol
        clear vert_exact_transverse
        u_sol = NaN(1,Mesh(ni).n_vert);
        %t_sol = NaN(1,Surface_vis(ni).n_vert);
        disp_probs = false;
        parameters_transverse = NaN(length(u_sol),2);
        dist = u_sol; v_sol = dist; 
      
        transverse_vel = NaN(Mesh(ni).n_vert,3);
        max_transverse_tries = 50;
        %         dist_tol = 3E-1;  % 1E-1 now seems to small?
        dist_tols = [    0.1  0.5 1 2  ];
        tries = zeros(size(dist));
        parameter_tol = 1E-2;
        %     tol = 0.1;
        param_err = NaN(Mesh(ni).n_vert,2);
        transverse_dist_error = NaN(length(u_sol),1);
        transverse_param_error = NaN(length(u_sol),2);
        disp('starting transverse');
        %         meshbad = false;
        if progress_monitors
            ppm = ParforProgMon('Transverse submesh    ', length(u_sol));
        else
            ppm = [];
        end
        
        parfor (vert_i = 1:length(u_sol), input.performance.nthreads)
            %                                     for vert_i = 1:length(u_sol)
            
            
            
            %               disp(['transverse vert ',num2str(vert_i),'     ','doneness = ',num2str(vert_i / Mesh(ni).n_vert)]);
            %  vert_i / Mesh(ni).n_vert
            vert = Mesh(ni).verts(vert_i,:);
            %   vert = Surface_vis(ni).verts(vert_i,:);
            
            distance = @(x) sqrt(sum((vert' - transverse_parameterized(x(1),x(2), time, geom.transverse)).^2));
            % answer = fmincon(distance,[0,0],[],[],[],[],[0 0],[geom.transverse.u_max geom.transverse.v_max],[],optimset('maxfunevals',1E20,'maxiter',1E20,'tolx',1E-13,'tolfun',1E-13)); %limit upper bound by t_max so points on hemisphere will fall at this parameter value
            dist(vert_i) = Inf; u_sol(vert_i) = Inf;  v_sol(vert_i) = Inf;
            
            if static_transverse
                still_bad = false;
            else
                still_bad = true;
            end
            
            for d = 1:length(dist_tols)
                if ~still_bad  %success
                    break
                end
                
                dist_tol = dist_tols(d);
                
                %                 if d == length(dist_tols)
                %                     disp('trying final dist tol')
                %                     %                         pause
                %                 end
                
                tc = 0;  still_bad = false;
                
                
                while dist(vert_i) > dist_tol || u_sol(vert_i) < 0 - parameter_tol || v_sol(vert_i) < 0 - parameter_tol || u_sol(vert_i) > geom.transverse.u_max + parameter_tol || v_sol(vert_i) > geom.transverse.v_max + parameter_tol
                    %                 answer = fminsearch(distance,[geom.transverse.u_max*rand, geom.transverse.v_max*rand],optimset('maxfunevals',1E2,'maxiter',1E3,'tolx',1E-6,'tolfun',1E-13))
                    
                    tc = tc + 1;
                    
                    if tc > max_transverse_tries
                        still_bad = true;
                        
                        break
                    end
                    
                    
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
                    
                    %                     if tc > max_transverse_tries
                    %                         meshbad = meshbad | true;
                    %                         error('Max transverse tries reached');
                    %                     end
                    
                end  %while
                
            end %dist_tols
            
            if d == length(dist_tols) && still_bad
                %                 meshbad = true;
                error(['Problem with transverse submesh vert_i ',num2str(vert_i)]);
            end
            if progress_monitors
                
                ppm.increment();
            end
        end  %parfor verts
        
        if progress_monitors && exist('ppm','var')
            
            delete(ppm);
        end
        %         if meshbad
        %             error('Problem with transverse submesh');
        %         end
        
        transverse_dist_error(f) = max(dist);
        transverse_param_error(f,:) = max(param_err);
        
        %         if transverse_dist_error(f) > 0.1
        %             disp('stopafra')
        %             pause
        %         end
        timings.transverse = toc;
        
        if ~static_transverse
        Mesh(ni).verts = vert_exact_transverse;
        else
            transverse_vel(:) = 0;
        end
        
        for nii = setdiff(1:length(Mesh),ni)
            [is_shared, inds] = ismember( Mesh(nii).indices.orig.vert , Mesh(ni).indices.orig.vert );
            
            Mesh(nii).verts( is_shared  , :) = Mesh(ni).verts( inds(is_shared), :);
        end
    end
    
    
    %%
    clear hair_vels hair_parameters
    hair_cases = {'Coplanar_Hairs','Normal_Top_Hairs','Normal_Bottom_Hairs'};
    for hc = 1:length(hair_cases)
        hair_case = hair_cases{hc};
        
        if ismember(hair_case,names)
            tic;
            temp = {Mesh.name};
            ni = find(strcmp(temp,hair_case));
            clear dist u_sol h_sol
            u_sol = NaN(1,Mesh(ni).n_vert);
            %t_sol = NaN(1,Surface_vis(ni).n_vert);
            disp_probs = false;
            
            dist = u_sol; h_sol = dist; wingtip_vel = NaN(Mesh(ni).n_vert,3);
            parameters_wingtip = NaN(length(u_sol),2);
            %dist_tol = 5E-2;
            %         dist_tol = 1E-1;  %sometimes get 0.4 off...
            
            %dist_tols = [0.05 0.1 0.15 0.2];  %parallel_2
            dist_tols = [    0.1    0.5 1 2 ];  %parallel_2
            %          dist_tol = 0.04;  %parallel_1
            %          dist_tol = 0.2;
            max_tries = 50;  % 5 is too small, sometimes fails
            tries = zeros(size(dist));
            parameter_tol = 5E-2;
            %     tol = 0.1;
            clear vert_exact_wingtip
            param_err = NaN(Mesh(ni).n_vert,2);
            wingtip_dist_error = NaN(length(u_sol),1);
            wingtip_param_error = NaN(length(u_sol),2);
            disp(['starting ',hair_case]);
            %
            if ~isempty(last_solution.(hair_case))
                NS = createns(last_solution.(hair_case).verts);
            end
            
            
            meshbad = false;
            if progress_monitors
                ppm = ParforProgMon([hair_case,'    '], length(u_sol));
            else
                ppm = [];
            end
            parfor (vert_i = 1:length(u_sol), input.performance.nthreads)
                %                             for vert_i = 1:length(u_sol)
                answer = [];
                guess0 = [];
                
                
                %   for vert_i = 1:length(u_sol)
                %                   disp([hair_case,' vert ',num2str(vert_i),'     ','doneness ',num2str( vert_i / Mesh(ni).n_vert) ]);
                %  vert_i / Mesh(ni).n_vert
                vert = Mesh(ni).verts(vert_i,:);
                %   vert = Surface_vis(ni).verts(vert_i,:);
                if  ~isempty(last_solution.(hair_case))
                    Idx = knnsearch(NS,vert);
                    guess0 = last_solution.(hair_case).parameters(Idx,:);
                end
                
                
                distance = @(x) sqrt(sum((vert' - transverse_hairs_parameterized(x(1),x(2), time, geom.transverse , hair_case)).^2));
                % answer = fmincon(distance,[0,0],[],[],[],[],[0 0],[geom.transverse.u_max geom.transverse.v_max],[],optimset('maxfunevals',1E20,'maxiter',1E20,'tolx',1E-13,'tolfun',1E-13)); %limit upper bound by t_max so points on hemisphere will fall at this parameter value
                dist(vert_i) = Inf; u_sol(vert_i) = Inf;  h_sol(vert_i) = Inf;
                
                if static_transverse
                    still_bad = false;
                else
                still_bad = true;
                end
                
                
                for d = 1:length(dist_tols)
                    if ~still_bad  %success
                        break
                    end
                    
                    dist_tol = dist_tols(d);
                    
                    if d == length(dist_tols)
                        %                         disp('trying final dist tol')
                        %                         pause
                    end
                    
                    tc = 0;  still_bad = false;
                    while dist(vert_i) > dist_tol || u_sol(vert_i) < 0 - parameter_tol || h_sol(vert_i) < geom.hairs.(hair_case).h_min - parameter_tol || u_sol(vert_i) > geom.transverse.u_max + parameter_tol || h_sol(vert_i) > geom.hairs.(hair_case).h_max + parameter_tol
                        % answer = fminsearch(distance,[geom.transverse.u_max*rand, geom.transverse.h_min + (geom.transverse.h_max - geom.transverse.h_min)*rand],optimset('maxfunevals',1E20,'maxiter',1E4,'tolx',1E-13,'tolfun',1E-13));
                        tc = tc + 1;
                        
                        if tc > max_tries
                            still_bad = true;
                            
                            break
                        end
                        
                        if tc == 1 && ~isempty(last_solution.(hair_case))
                            guess = guess0;
                        else
                            %                             disp('guess0 didn''t work')
                            guess =  [geom.transverse.u_max*rand, geom.hairs.(hair_case).h_min + (geom.hairs.(hair_case).h_max - geom.hairs.(hair_case).h_min)*rand];
                        end
                        
                        %                 answer = fmincon(distance, [geom.transverse.u_max*rand, geom.transverse.h_min + (geom.transverse.h_max - geom.transverse.h_min)*rand], [],[],[],[],[0 geom.transverse.h_min]',[geom.transverse.u_max geom.transverse.h_max]',[],optimoptions('fmincon','display','off'))
                        %                         if meshbad
                        %                             answer = NaN(1,2);
                        %                         else
                        try
                            answer = patternsearch(distance,guess, [],[],[],[],[0 geom.hairs.(hair_case).h_min]',[geom.transverse.u_max geom.hairs.(hair_case).h_max]',[],optimoptions('patternsearch','display','off','PollOrderAlgorithm', 'random', 'UseCompletePoll', false,'AccelerateMesh',true,'cache','on'));
                            %                         answer = fmincon(      distance, [geom.transverse.u_max*rand, geom.transverse.h_min + (geom.transverse.h_max - geom.transverse.h_min)*rand], [],[],[],[],[0 geom.transverse.h_min]',[geom.transverse.u_max geom.transverse.h_max]',[],optimoptions('fmincon','display','off'));
                            % answer
                        catch
                            answer = NaN(1,2);
                        end
                        %                         end
                        
                        temp =  transverse_hairs_parameterized(answer(1),answer(2), time, geom.transverse , hair_case);
                        
                        vert_exact_wingtip(vert_i,:) = temp;
                        
                        if any(isnan(  temp ) )
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
                        
                        if h_sol(vert_i) < geom.hairs.(hair_case).h_min
                            h_err = geom.hairs.(hair_case).h_min - h_sol(vert_i);
                        elseif h_sol(vert_i) > geom.hairs.(hair_case).h_max
                            h_err = geom.hairs.(hair_case).h_max - h_sol(vert_i);
                        else
                            h_err = 0;
                        end
                        
                        param_err(vert_i,:) = [u_err h_err];
                        
                        [~,wingtip_vel(vert_i,:)] = transverse_hairs_parameterized(u_sol(vert_i),h_sol(vert_i), time, geom.transverse , hair_case);
                        tries(vert_i) =  tries(vert_i) + 1;
                        %   [u_sol(vert_i)  v_sol(vert_i)]
                        % dist(vert_i)
                        if   dist(vert_i) > dist_tol && disp_probs
                            disp(['vert ',num2str(vert_i),'     dist > dist_tol, dist = ',num2str(dist(vert_i))]);
                        end
                        
                        if u_sol(vert_i) < 0 - parameter_tol  && disp_probs
                            disp(['u_sol < 0, u_sol = ',num2str(u_sol(vert_i))]);
                        end
                        
                        if h_sol(vert_i) < geom.hairs.(hair_case).h_min - parameter_tol  && disp_probs
                            disp(['h_sol < ',num2str(geom.hairs.(hair_case).h_min),', h_sol = ',num2str(h_sol(vert_i))]);
                        end
                        
                        if u_sol(vert_i) > geom.transverse.u_max + parameter_tol   && disp_probs
                            disp(['u_sol > u_max, u_sol - u_max = ',num2str(u_sol(vert_i) - geom.transverse.u_max)]);
                        end
                        
                        if h_sol(vert_i) > geom.hairs.(hair_case).h_max + parameter_tol  && disp_probs
                            disp(['h_sol > ',num2str(geom.hairs.(hair_case).h_max),', h_sol = ',num2str(h_sol(vert_i))]);
                        end
                        
                        
                        
                    end  %while
                    
                end  %for dist_tols
                
                if d == length(dist_tols) && still_bad
                    meshbad = meshbad | true;
                    %                     disp(['vert ',num2str(vert_i),' is bad']);
                    error(['Problem with ',hair_case,' vert_i ',num2str(vert_i)]);
                end
                if progress_monitors
                    ppm.increment();
                    
                end
            end  %parfor verts
            if progress_monitors && exist('ppm','var')
                delete(ppm);
            end
            %             if meshbad
            %                 %                 pause
            %                 error('Problem with hair submesh');
            %             end
            %%
            wingtip_dist_error(f) = max(dist);
            wingtip_param_error(f,:) = max(param_err);
            
            %         if wingtip_dist_error(f) > dist_tol
            %             disp('stopafra')
            %             pause
            %         end
            timings.wingtip = toc;
            if ~static_transverse
            Mesh(ni).verts = vert_exact_wingtip;
            else
                wingtip_vel(:) = 0;
            end
            
            for nii = setdiff(1:length(Mesh),ni)
                [is_shared, inds] = ismember( Mesh(nii).indices.orig.vert , Mesh(ni).indices.orig.vert );
                
                Mesh(nii).verts( is_shared  , :) = Mesh(ni).verts( inds(is_shared), :);
            end
            
            hair_vels.(hair_case) = wingtip_vel;
            hair_parameters.(hair_case) = parameters_wingtip;
            if static_transverse
                last_solution.(hair_case).verts = Mesh(ni).verts;
            else
            last_solution.(hair_case).verts = vert_exact_wingtip;
            end
            last_solution.(hair_case).parameters = parameters_wingtip;
            
        end  % if we have this hair case
        
        
    end  % hair cases loop (coplanar, normal top, normal bottom)
    %%
    if meshbad
        continue
    end
    
    all_verts = vertcat(Mesh.verts);
    if any(isnan(all_verts(:)))
        stoperdoodle
    end
    
    %         stoperdoodle
    %%
    clear Metadata
    Metadata.geom = geom;
    Metadata.time = time;
    
    for i = 1:length(names)  %ensures order within Metadata matches order in which submeshes were loaded (which is order of names_L)
        name = names{i};
        switch name
            case 'Body'
                Metadata.Body.BCs = zeros(Mesh(i).n_vert,3);
                
            case 'Transverse'
%                 if static_transverse
%                 Metadata.Transverse.BCs = zeros(Mesh(i).n_vert,3);
%                 Metadata.Transverse.parameters = [];
%                 else
                    Metadata.Transverse.BCs = transverse_vel;
                Metadata.Transverse.parameters = parameters_transverse;
%                 end
            case 'Tail'
                Metadata.Tail.BCs = tail_vel;
                Metadata.Tail.parameters = parameters_tail;
            case 'Coplanar_Hairs'
%                 if static_transverse
%                 Metadata.Coplanar_Hairs.BCs = zeros(Mesh(i).n_vert,3);
%                 Metadata.Coplanar_Hairs.parameters =[];
%                 else
                    Metadata.Coplanar_Hairs.BCs = hair_vels.Coplanar_Hairs;
                Metadata.Coplanar_Hairs.parameters = hair_parameters.Coplanar_Hairs;
%                 end
            case 'Normal_Top_Hairs'
%                 if static_transverse
%                 Metadata.Normal_Top_Hairs.BCs = zeros(Mesh(i).n_vert,3);
%                 Metadata.Normal_Top_Hairs.parameters = [];
%                 else
                      Metadata.Normal_Top_Hairs.BCs = hair_vels.Normal_Top_Hairs;
                Metadata.Normal_Top_Hairs.parameters = hair_parameters.Normal_Top_Hairs;
%                 end
            case 'Normal_Bottom_Hairs'
%                 if static_transverse
%                 Metadata.Normal_Bottom_Hairs.BCs =  zeros(Mesh(i).n_vert,3);
%                 Metadata.Normal_Bottom_Hairs.parameters = [];
%                 else
                     Metadata.Normal_Bottom_Hairs.BCs = hair_vels.Normal_Bottom_Hairs;
                Metadata.Normal_Bottom_Hairs.parameters = hair_parameters.Normal_Bottom_Hairs; 
%                 end
        end
        Metadata.(name).indices.orig.vert = Mesh(i).indices.orig.vert;
        %         Metadata.(name).rand_inds = Metadata_temp(i).mesh.rand_inds;
        %  Metadata.(name).rand_inds = [1:Mesh(i).n_vert];
    end
    
    ind1 = strfind(Files{1}{f},'step_');
    save([folder_processed,'Metadata_',Files{1}{f}(ind1:end-4),'.mat'],'Metadata');
    
    
    for i = 1:length(names)
        
        
        filename = [folder_processed, names{i},'_',Files{1}{f}(ind1:end)];
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