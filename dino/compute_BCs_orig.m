folder = 'E:\Hull\dinoflagellate\meshes_stub_groups\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_stub_groups\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_grouped_fast\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail_nosphere2\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail_nosphere_tail_angle_0\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_thick_tail_nosphere_tail_angle_0\';
folder = 'E:\Hull\dinoflagellate\meshes_biggermin\';

input.performance.nthreads = 8;
store_constants = true;

clear Inds Files
names = {'Body','Transverse','Tail'};
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

if ~isequal(Times{1},Times{2},Times{3})
    disp('problemo')
    pause
else
    Time = Times{1};
end
%%
for f = 1:length(Files{1})
    
    disp(['Step ',num2str(f),' out of ',num2str(length(Files{1}))]);
    
    time = Time(f);
    clear Mesh rand_inds
   
    [Mesh(1),rand_inds(1)] = load_mesh([folder,Files{1}{f}],[],[],'mesh');  %save rand vertex permutation for use later, so all future calculations use same exact mesh!
    [Mesh(2),rand_inds(2)] = load_mesh([folder,Files{2}{f}],[],[],'mesh');
    [Mesh(3),rand_inds(3)] = load_mesh([folder,Files{3}{f}],[],[],'mesh');
   
   
    Mesh(1).name = 'body';
    Mesh(2).name = 'transverse';
    Mesh(3).name = 'tail';
    
    for si = 1:length(Mesh)
        Mesh(si).orientation = [1 0 0; 0 1 0; 0 0 1;]';
        Mesh(si).refpoints = [0 0 0]';  
    end
    
    [Mesh] = global_inds(Mesh);
    [Mesh] = renumber_Mesh(Mesh);
    
    if store_constants
        
        temp_input.performance.nthreads = 8;
        
        temp_input.accuracy.integration_tol.area.abstol = 0;
        temp_input.accuracy.integration_tol.area.reltol = 1000;
        temp_input.accuracy.integration_tol.area.maxevals = Inf;
        
        temp_input.accuracy.integration_tol.centroid.abstol = 0;
        temp_input.accuracy.integration_tol.centroid.reltol = 1000;
        temp_input.accuracy.integration_tol.centroid.maxevals = Inf;
        
        temp_input.accuracy.integration_tol.volume.abstol = 0;
        temp_input.accuracy.integration_tol.volume.reltol = 1000;
        temp_input.accuracy.integration_tol.volume.maxevals = Inf;
        
        [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
    end
    
    clear geom
    
    geom.lambda = 23.7;
    geom.nlambda = 1.5;
    geom.amp = [1.752  5.9];
    geom.kE = [3 * 0.24 0.1];
    geom.omega = 2 * pi / 0.0219 ;
    
    geom.t_transition = 6.01151778049;
    geom.t_transition = 5.91727000088;
     geom.t_transition = 2.51113977498;
      geom.t_transition =  2.51113977498;
     geom.t_transition =  3.29161918769;  %thin tail angled
     
    geom.t_transition =  3.20735258036;  %thick tail angled
    
     geom.t_transition = 2.85754556161;
%      geom.t_transition = 2.85754556335;
    
    geom.translation = ([1.31837142643 2.77555756156e-16 0.925] -   [0.0 -0.0 0.0])';
    geom.translation = ([1.30337142643 2.77555756156e-16 0.925] -   [0.0 -0.0 0.0])';
   
    geom.translation = ([ 6.04462993646 3.5527136788e-15 15.2] -   [0.0 -0.0 0.0])';
     geom.translation = ([ 6.04462993646 3.5527136788e-15 15.2])';
     
       geom.translation = ([ 6.36348374047 3.5527136788e-15 15.2])';  %thin tail angled
      geom.translation = ([ 6.04563244435 3.5527136788e-15 15.2])';   %thick tail angled
     
       geom.translation = [6.07001857152 4.60467196479e-15 14.0]';
%        geom.translation = [6.0700185781  4.60467196479e-15 14.0]';
   
    
    geom.rotation = [0 0 pi]';
    geom.radius = 0.15   * 1;
    geom.t_max = 2*pi * geom.nlambda;  %t value where hemispherical end cap is centered
    
    safety_factor = 1.5;  % t value at end of sphere must be less than safety_factor * geom.t_max
    
    geom.tail_angle = -25 *0  ;
    
    geom.rotation_pt = [6.85279653889 3.5527136788e-15 15.2]';  %thin tail angled
     geom.rotation_pt = [6.85279653889 3.5527136788e-15 15.2]';  %thick tail angled
     
       geom.rotation_pt  = [6.07001857152 4.60467196479e-15 14.0]';
     
    geom.rotation_vec = [0 0 1];  %always vertically up along z dir
    
    %time = 6.28318530718;
    
    names = {Mesh.name};
    ni = find(strcmp(names,'tail'));
    
    t_sol = NaN(1,Mesh(ni).n_vert);
    dist = t_sol; theta = dist;  tail_vel = NaN(Mesh(ni).n_vert,3);  n_t = tail_vel;
    pts = tail_vel;  V_surfaces = pts;  V_us = pts;
    
    
    for vert_i = 1:length(t_sol)
         vert_i / Mesh(ni).n_vert
        vert = Mesh(ni).verts(vert_i,:);
        %   vert = Surface_vis(ni).verts(vert_i,:);
        
        distance = @(t) sqrt(sum((vert - tail_parameterized(t, time, geom)).^2));  %distance between vert and tail centerline as function of t
        
        t_sol(vert_i) = fminbnd(distance,0,geom.t_max*safety_factor,optimset('maxfunevals',1E8,'maxiter',1E8,'tolx',1E-12)); %limit upper bound by t_max so points on hemisphere will fall at this parameter value
        
        dist(vert_i) = distance(t_sol(vert_i));
        
        %for spherical end cap, add vel of centerline pt to time der of vector
        %formed by rotating V_u in 3D around a new vector normal to both V_u and
        %V_surface (gotten with a cross product of the two) by an angle from the
        %dot product of the two
        
        [pt, vel, der, V_u] = tail_parameterized(t_sol(vert_i), time, geom);
        pts(vert_i,:) = pt;
        V_surface = vert - pt;
        
        V_surface = V_surface / sqrt(sum(V_surface.^2));
        V_u = V_u / sqrt(sum(V_u.^2));
        
        V_surfaces(vert_i,:) = V_surface;
        V_us(vert_i,:) = V_u;
        
        if vert(3) == 0
            stopafra
        end
        
        theta(vert_i) = acos(dot(V_surface,V_u)) * (- sign(V_surface(3)));
        V_rot = cross(V_u, V_surface);
        
        [n_t(vert_i,:)] = surface_vector_vel(t_sol(vert_i), theta(vert_i), time, geom);
        
        tail_vel(vert_i,:) = vel + dist(vert_i) * n_t(vert_i,:);
        
        if sqrt(sum(vel.^2)) < 1E-10
            %  n_t(vert_i,:)
            % pause
        end
        %V_surface
        % V_u
        % V(i,:) = [0 1 0; -1 0 0; 0 0 1] * der(i,:)';  %rotate tangent vector 90 deg around z-axis (which is normal to movement plane) to obtain "upward" vector
        
        % [pt, vel, der] = change_ref_frame(pt, vel, der, geom)
    end
    factor  = 0.99;  %to make sure we don't look at points on flat end surface at t = tmax
    
    tail_error = abs([min(dist(t_sol < geom.t_max*factor)) max(dist(t_sol < geom.t_max*factor))] - geom.radius)  %an estimate of error = difference between distances and actual radius
    
    if max(tail_error) > 1E-2
        disp('stopafra')
        pause
    end
    %%
%     t  = linspace(0,2*pi*geom.nlambda*1.01,3000);
%     inc = 12;
%     for Time = linspace(0,20,20)
%         [pt, vel, der] = tail_parameterized(t, Time, geom);
%     
%         figure(453)
%         plot3(pt(:,1),pt(:,2),pt(:,3),'-','markerfacecolor','k');
%         hold on
%         %quiver3(pt(1:inc:end,1),pt(1:inc:end,2),pt(1:inc:end,3),vel(1:inc:end,1),vel(1:inc:end,2),vel(1:inc:end,3), 0.1);
%         % q = quiver3(pt(1:inc:end,1),pt(1:inc:end,2),pt(1:inc:end,3),der(1:inc:end,1),der(1:inc:end,2),der(1:inc:end,3), 1);
%         % set(q,'showarrowhead','off');
%         [pt, vel, der] = tail_parameterized(geom.t_max, Time, geom);
%         plot3(pt(:,1),pt(:,2),pt(:,3),'o','markerfacecolor','k');
%     
%         hold off
%         title(num2str(Time));
%         grid on
%         %ylim([-0.5 0.5]);
%         view([0 90]);
%         pause
%     
%         drawnow
%     
%     end
    %%
%     clear Pts
%     tvec = linspace(0,geom.t_max,300);
%     for i = 1:length(tvec)
%      [pt, vel, der] = tail_parameterized(tvec(i), 0, geom);
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

    
    %%
    
    geom.R = 14.8;
    geom.R = 13.05;
    
    %geom.R = 15.7;  big groove I think
    geom.d = 0;
    %geom.w = 0.1;
    geom.w = 0.6;  %4
    geom.rf = 1.8;
    geom.c = 14;
    geom.omega = 2 * pi / 0.0258;
    
    geom.b = 3.4 ;  %wavelength of underlying helix - 0 for no helix, just a circle
    geom.N_revs = 1;
%     geom.u_max = 2*pi*geom.N_revs;
%     geom.u_max = 6.151446974906251;  %no longer based on N_revs since sheet stops where groove ends, which Salome determines
    geom.u_max = 6.15144695058772;
    geom.u_max = 6.15144697491;
    
    geom.v_max = 1;
%     geom.shift = 8.7;
    geom.shift = 7.5*0.825*1.25;
    
    
    stopooo
    
    
    %% generate theoretical points more or less evenly spaced over transverse
    
    u = linspace(0,geom.u_max,5000);
    v = linspace(0,geom.v_max,100);
    c = 0;
    pts = NaN(length(u)*length(v),3);  U = NaN(length(u)*length(v),1); V = U;
    for i = 1:length(u)
        i/length(u)
    
        for j = 1:length(v)
            c = c+1;
            U(c) = u(i);  V(c) = v(j);
    
            pts(c,:) = transverse_parameterized(u(i),v(j), time, geom)';
        end
    end
    
    
    
    %% show closest points on theoretical surface to each vertex of actual mesh
    u = u_sol;  v = v_sol;  %u and v had better be same size since they match with verts
    
       %c = 0;
    pts_fit = NaN(length(u),3);  U = NaN(length(u),1); V = U;
    parfor i = 1:length(u)
      %  i/length(u)

           % c = c+1;
%             U(c) = u(i);  V(c) = v(j);
    
            pts_fit(i,:) = transverse_parameterized(u(i),v(i), time, geom)';
      
    end
    %%
    clear th
    for i = 1:size(pts_fit,1)
        th(i) = text(pts_fit(i,1),pts_fit(i,2),pts_fit(i,3),num2str(v_sol(i)),'fontsize',12);
    end
    %%
    
    names = {Mesh.name};
    ni = find(strcmp(names,'transverse'));
    clear dist u_sol v_sol
    u_sol = NaN(1,Mesh(ni).n_vert);
    %t_sol = NaN(1,Surface_vis(ni).n_vert);
    disp_probs = true;
    
    dist = u_sol; v_sol = dist; transverse_vel = NaN(Mesh(ni).n_vert,3);
    dist_tol = 5E-2;
%     dist_tol = 0.1;
    tries = zeros(size(dist));
    tol = 1E-2;
%     tol = 0.1;
    
    parfor vert_i = 1:length(u_sol)
        vert_i / Mesh(ni).n_vert
        vert = Mesh(ni).verts(vert_i,:);
        %   vert = Surface_vis(ni).verts(vert_i,:);
        
        distance = @(x) sqrt(sum((vert' - transverse_parameterized(x(1),x(2), time, geom)).^2));  
        % answer = fmincon(distance,[0,0],[],[],[],[],[0 0],[geom.u_max geom.v_max],[],optimset('maxfunevals',1E20,'maxiter',1E20,'tolx',1E-13,'tolfun',1E-13)); %limit upper bound by t_max so points on hemisphere will fall at this parameter value
        dist(vert_i) = Inf; u_sol(vert_i) = Inf;  v_sol(vert_i) = Inf;
        
        while dist(vert_i) > dist_tol || u_sol(vert_i) < 0 - tol || v_sol(vert_i) < 0 - tol || u_sol(vert_i) > geom.u_max + tol || v_sol(vert_i) > geom.v_max + tol
            answer = fminsearch(distance,[geom.u_max*rand, geom.v_max*rand],optimset('maxfunevals',1E20,'maxiter',1E20,'tolx',1E-13,'tolfun',1E-13));
            
            u_sol(vert_i) = answer(1);  v_sol(vert_i) = answer(2);
            
            dist(vert_i) = distance(answer);
            
            [~,transverse_vel(vert_i,:)] = transverse_parameterized(u_sol(vert_i),v_sol(vert_i), time, geom);
            tries(vert_i) =  tries(vert_i) + 1;
          %   [u_sol(vert_i)  v_sol(vert_i)]
           % dist(vert_i)
          if   dist(vert_i) > dist_tol && disp_probs
              disp(['dist > dist_tol, dist = ',num2str(dist(vert_i))]);
          end
          
          if u_sol(vert_i) < 0 - tol  && disp_probs
              disp(['u_sol < 0, u_sol = ',num2str(u_sol(vert_i))]);
          end
          
          if v_sol(vert_i) < 0 - tol  && disp_probs
                disp(['v_sol < 0, v_sol = ',num2str(v_sol(vert_i))]);
          end
              
          if u_sol(vert_i) > geom.u_max + tol   && disp_probs
                 disp(['u_sol > u_max, u_sol - u_max = ',num2str(u_sol(vert_i) - geom.u_max)]);
          end
          
          if v_sol(vert_i) > geom.v_max + tol  && disp_probs
                 disp(['v_sol > v_max, v_sol - v_max = ',num2str(v_sol(vert_i) - geom.v_max)]);
          end
          
            
            
        end
        
        
        
    end
    
    transverse_error = max(dist)
    if transverse_error > 1E-2
        disp('stopafra')
        pause
    end
    %%
    clear Metadata
    Metadata.geom = geom;
    Metadata.time = time;
    Metadata.body.BCs = zeros(Mesh(1).n_vert,3);
    Metadata.transverse.BCs = transverse_vel;
    Metadata.tail.BCs = tail_vel;
    Metadata.body.indices.orig.vert = Mesh(1).indices.orig.vert;
    Metadata.transverse.indices.orig.vert = Mesh(2).indices.orig.vert;
    Metadata.tail.indices.orig.vert = Mesh(3).indices.orig.vert;
    Metadata.body.rand_inds = rand_inds(1).mesh.rand_inds;
    Metadata.transverse.rand_inds = rand_inds(2).mesh.rand_inds;
    Metadata.tail.rand_inds = rand_inds(3).mesh.rand_inds;
    

    save([folder,'Metadata_',Files{1}{f}(6:end-4),'.mat'],'Metadata');
    
 stopp
    
end  %file loop

stopa


%%

try, delete(q1), end;   try, delete(q2), end;
inc = 1;  ni = 3;
q1 = quiver3(Mesh(ni).verts(1:inc:end,1),Mesh(ni).verts(1:inc:end,2),Mesh(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3), 1);
%q = quiver3(Surface_vis(ni).verts(1:inc:end,1),Surface_vis(ni).verts(1:inc:end,2),Surface_vis(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3));
set(q1,'color','k');


inc = 1;  ni = 2;
q2 = quiver3(Mesh(ni).verts(1:inc:end,1),Mesh(ni).verts(1:inc:end,2),Mesh(ni).verts(1:inc:end,3),transverse_vel(1:inc:end,1),transverse_vel(1:inc:end,2),transverse_vel(1:inc:end,3), 1);
%q = quiver3(Surface_vis(ni).verts(1:inc:end,1),Surface_vis(ni).verts(1:inc:end,2),Surface_vis(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3));
set(q2,'color','k');


%%

inc = 1;  ni = 3;
q1 = quiver3(Mesh(ni).verts(1:inc:end,1),Mesh(ni).verts(1:inc:end,2),Mesh(ni).verts(1:inc:end,3),BCs.freeswim.tail(1:inc:end,1),BCs.freeswim.tail(1:inc:end,2),BCs.freeswim.tail(1:inc:end,3), 1);
%q = quiver3(Surface_vis(ni).verts(1:inc:end,1),Surface_vis(ni).verts(1:inc:end,2),Surface_vis(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3));
set(q1,'color','k');


inc = 1;  ni = 2;
q2 = quiver3(Mesh(ni).verts(1:inc:end,1),Mesh(ni).verts(1:inc:end,2),Mesh(ni).verts(1:inc:end,3),BCs.freeswim.transverse(1:inc:end,1),BCs.freeswim.transverse(1:inc:end,2),BCs.freeswim.transverse(1:inc:end,3), 0.5);
%q = quiver3(Surface_vis(ni).verts(1:inc:end,1),Surface_vis(ni).verts(1:inc:end,2),Surface_vis(ni).verts(1:inc:end,3),tail_vel(1:inc:end,1),tail_vel(1:inc:end,2),tail_vel(1:inc:end,3));
set(q2,'color','k');

