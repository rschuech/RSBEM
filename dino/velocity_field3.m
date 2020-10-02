
field_input.accuracy.integration_tol.field_vel.abstol = 1E-7;
field_input.accuracy.integration_tol.field_vel.abstol = 1E-6;
%field_input.accuracy.integration_tol.field_vel.abstol = 1E-5;

field_input.accuracy.integration_tol.field_vel.reltol = 0;
field_input.accuracy.integration_tol.field_vel.maxevals = Inf;
field_input.accuracy.eps2 = input.accuracy.eps2;
field_input.constants.multfactor = input.constants.multfactor;
field_input.constants.mu = input.constants.mu;

%integrand_constants = struct('eps2',assembly_input.accuracy.eps2,'x_col',Mesh(i_mesh_vert).verts(local_vert,:)');
% z = near(0,Z2);
% x = near(9,X2); %middle of groove  % 9
% y = near(0,Y2);
% load('C:\Users\rudi\Desktop\RD\dino_dumps_hairs\dino_parallel2_3_all_tailed3_body-transverse-tail-wingtip_dump.mat');
% folder  = 'C:\Users\rudi\Desktop\RD\dino dumps potato\';

% folder = 'C:\Users\rudi\Desktop\RD\dino_dumps_hairs_2_1.5_phases\';
% folder = 'C:\Users\rudi\Desktop\RD\dino_dumps_hairs_2_1.5_no_offset\';
folder = 'C:\Users\rudi\Desktop\RD\dino_dumps_hairs_3.5_1.5_phases\';

dump_folder = 'C:\Users\rudi\Desktop\RD\base_case_2_1.5\';
mesh_folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\regular\final\';

outfolder = 'C:\Users\rudi\Desktop\RD\dino velocity fields\hairs 2 1.5\';

do_rotation = true;
rotation_angle = [-47.5*pi/180  0  0]; % Z Y X  for z slice
% rotation_angle = [180*pi/180  0  0]; % Z Y X  for y slice
%  rotation_angle = [0  0  0]; % Z Y X  for x slice



dumps = dir([dump_folder,'*dump.mat']);
dumps = {dumps.name};

for d = 1 %1:length(dumps)
    [d length(dumps)];
    dump = dumps{d}
    
    
    load([dump_folder,dump],'Mesh_files','Solutions','input');
    if exist('mesh_folder','var') && ~isempty(mesh_folder)
        input.paths.datfolder = mesh_folder;
    end
    
    clear vel_fields
    % load('C:\Users\rudi\Desktop\RD\dino_dumps_hairs\dino_parallel2_3_all_tailed3_wingtip_dump.mat');
    % load('C:\Users\rudi\Desktop\RD\dino_dumps_hairs\dino_parallel2_3_all_tailed3_transverse_dump.mat');
    % load('C:\Users\rudi\Desktop\RD\dino_dumps_hairs\dino_parallel2_3_all_tailed3_body-transverse-tail-wingtip_dump.mat');
    
    
    % has Solutions struct
    Solutions_fields = setdiff( fieldnames(Mesh_files) , {'phase','time','metadata'});  %order of submeshes as far as Solutions (e.g. traction f) is concerned
    % hmm, it actually appears that the order of Solutions_fields is the same
    % as the order of Solutions.rand_inds, so probably don't need to do
    % anything to match the orders
    
    %%
    Xlist = [];  Ylist = [];  Zlist = [];
    
    spacing = 2.5;  %put vectors every this many microns  Lasse's images appear to have 3-4 vectors / 10 um
    
    Lx = [-75 70];
    Ly = [-65 65];
    Lz = [-65 65];
    
    Lx = [-105 120];
    Ly = [-100 65];
    Lz = [-60 60];
    
    
    
    
    x = 30.05; %top of groove % 12.5
    y = Ly(1):spacing:Ly(2);
    z = Lz(1):spacing:Lz(2);
    [X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
    Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];
    
    
    % horizontal planes parallel to tail beating plane
    
    z = 0;  %middle of body, below tail % 0
    x = Lx(1):spacing:Lx(2);
    y = Ly(1):spacing:Ly(2);
    [X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
    Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];
    
    
    % vertical, along tail axis, down through body
    
    y = 0;  %through center of body, vertically, through tail % 0
    x = Lx(1):spacing:Lx(2);
    z = Lz(1):spacing:Lz(2);
    [X,Y,Z] = meshgrid(x,y,z);  %meshgrid BS
    Xlist = [Xlist; X(:)];  Ylist = [Ylist; Y(:)];   Zlist = [Zlist; Z(:)];
    
    
    points = [Xlist(:) Ylist(:) Zlist(:)];
    points = unique(points,'rows');
    
    %         if phase_ind == 1
    
    %         end
    
    x = sort(unique([30.05 Lx(1):spacing:Lx(2)]));
    y = sort(unique([0 Ly(1):spacing:Ly(2)]));
    z = sort(unique([0 Lz(1):spacing:Lz(2)]));
    [X,Y,Z] = meshgrid(x,y,z);
    %%
    Meshes = [];
    
    % average over a phase cycle
    fh = waitbar(0,'phase cycle velocity field avg');
    for phase_ind = 1:length(Solutions.phase)
        waitbar(phase_ind / length(Solutions.phase) );
        disp(['On phase ind ',num2str(phase_ind),' of ',num2str(length(Solutions.phase))]);
        % phase_ind = 1;
        phase = Solutions.phase(phase_ind);
        
        [Mesh, Metadata, matrix_props, Metadata_rand_inds] = load_dino_mesh(phase, Mesh_files, input, Solutions.rand_inds{phase_ind});
        
        if do_rotation
            [Mesh,rotmat] = rotateMesh(Mesh, rotation_angle);
            % Metadata apparently never used here so don't bother rotating
            % it
            %             for si = 1:length(Solutions.f)
            Solutions.f{phase_ind} = reshape( rotmat * reshape( Solutions.f{phase_ind} , 3, [])  , [] , 1);  % reshape tracton vector to a matrix, rotate it, reshape back to a vector
            %             end
            Solutions.U(phase_ind,:) = (rotmat * Solutions.U(phase_ind,:)')';
            Solutions.Omega(phase_ind,:) = (rotmat * Solutions.Omega(phase_ind,:)')';
        end
        
        Meshes{phase_ind} = Mesh;
    
        
        tic;
        [field_vel] = field_velocity_mexed(points, Solutions.f{phase_ind}, Mesh, field_input);
        disp(['field_velocity_mexed took ',num2str(toc)]);
        
        U = NaN(size(X));  V = U;  W = U;
        parfor i = 1:numel(X)
            
            point = [X(i) Y(i) Z(i)];  %point in full 3D grid
            ind = find(  points(:,1) == point(1) & points(:,2) == point(2) & points(:,3) == point(3) );
            if ~isempty(ind)
                if numel(ind) > 1
                    stopafra
                end
                
                U(i) = field_vel(ind,1);  V(i) = field_vel(ind,2);  W(i) = field_vel(ind,3);
                %    pause
            else
                U(i) = NaN;  V(i) = NaN;  W(i) = NaN;  %no idea why this is needed - parfor is forgetting the initialization of U,V,W and defaulting it to zeros!
            end
        end
        
        vel_fields(phase_ind).phase = Solutions.phase(phase_ind);
        vel_fields(phase_ind).X = X; vel_fields(phase_ind).Y = Y;  vel_fields(phase_ind).Z = Z;
        vel_fields(phase_ind).U = U; vel_fields(phase_ind).V = V;  vel_fields(phase_ind).W = W;
        vel_fields(phase_ind).field_vel = field_vel;
        vel_fields(phase_ind).points = points;
        
    end
    close(fh)
    %%
    fh = waitbar(0,'conversion to body frame');
    for phase_ind = 1:length(vel_fields)
        waitbar(phase_ind / length(Solutions.phase) );
        field_vel = vel_fields(phase_ind).field_vel;
        points = vel_fields(phase_ind).points;
        omega_cross_r = cross( repmat(Solutions.Omega(phase_ind,:),size(points,1),1)    ,   points - repmat(Mesh(1).refpoints(:,1)',size(points,1),1)   ) ;
        temp = field_vel - repmat(Solutions.U(phase_ind,:),size(field_vel,1),1) - omega_cross_r;
        
        U_bodyframe = NaN(size(X)); V_bodyframe = U_bodyframe; W_bodyframe = U_bodyframe;
        parfor i = 1:numel(X)
            
            point = [X(i) Y(i) Z(i)];  %point in full 3D grid
            ind = find(  points(:,1) == point(1) & points(:,2) == point(2) & points(:,3) == point(3) );
            if ~isempty(ind)
                if numel(ind) > 1
                    stopafra
                end
                
                U_bodyframe(i) = temp(ind,1); V_bodyframe(i) = temp(ind,2);  W_bodyframe(i) = temp(ind,3);
                %    pause
            else
                U_bodyframe(i) = NaN; V_bodyframe(i) = NaN;  W_bodyframe(i) = NaN;  %no idea why this is needed - parfor is forgetting the initialization of U,V,W and defaulting it to zeros!
            end
        end
        
        vel_fields(phase_ind).bodyframe.U = U_bodyframe;   vel_fields(phase_ind).bodyframe.V = V_bodyframe;   vel_fields(phase_ind).bodyframe.W = W_bodyframe;
    end
    close(fh)
    
    %%
    
%     interpolants_U = struct('x0',cell(size(vel_fields(1).X)));  interpolants_U_bodyframe = interpolants_U;
%     interpolants_V = struct('x0',cell(size(vel_fields(1).X)));  interpolants_V_bodyframe = interpolants_V;
%     interpolants_W = struct('x0',cell(size(vel_fields(1).X)));  interpolants_W_bodyframe = interpolants_W;
    names = {'U','V','W'};
    vel_field_avg_U = NaN(size(U));  vel_field_avg_V = NaN(size(V));  vel_field_avg_W = NaN(size(W));
    vel_field_avg_U_bodyframe = NaN(size(U));  vel_field_avg_V_bodyframe = NaN(size(V));  vel_field_avg_W_bodyframe = NaN(size(W));
    
    parfor i = 1:numel(vel_fields(1).X)
        temp = []; temp2 = []; fun = [];  fun2 = [];
        for n = 1:length(names)
            j = [];
            y = []; yb = [];
            for j = 1:length(Solutions.phase)
                y(j) = vel_fields(j).(names{n})(i);
                yb(j) = vel_fields(j).bodyframe.(names{n})(i);
            end
            switch input.kinematics_interp_method
                case 'trig'
                    temp = trig_interp_fit(Solutions.phase,y);
                    temp2 = trig_interp_fit(Solutions.phase,yb);
                    
%                     fields = fieldnames(temp);
%                     for f = 1:length(fields)
%                         switch names{n}
%                             case 'U'
%                                 [interpolants_U(i).(fields{f})] = temp.(fields{f});
%                                 [interpolants_U_bodyframe(i).(fields{f})] = temp2.(fields{f});
%                                 
%                             case 'V'
%                                 [interpolants_V(i).(fields{f})] = temp.(fields{f});
%                                 [interpolants_V_bodyframe(i).(fields{f})] = temp2.(fields{f});
%                             case 'W'
%                                 [interpolants_W(i).(fields{f})] = temp.(fields{f});
%                                 [interpolants_W_bodyframe(i).(fields{f})] = temp2.(fields{f});
%                         end
%                     end
                    
                case 'spline'
                    
                    filter = ~ismember(roundn(Solutions.phase,input.interp_phases_tol), roundn(input.discard_interp_phases,input.interp_phases_tol));
                   if ~all(isnan(y(filter)))
                    temp = csape([Solutions.phase(filter) ; phase_bounds(2)]' , [y(filter)  y(1)],'periodic'); % assume that first entry is phase = 0, which is never discarded
                    temp2 = csape([Solutions.phase(filter) ; phase_bounds(2)]' , [yb(filter)  yb(1)],'periodic');
                   else
                       temp = [];  temp2 = [];
                   end
%                     switch names{n}
%                         case 'U'
%                             [interpolants_U(i)] = temp;
%                             [interpolants_U_bodyframe(i)] = temp2;
%                             
%                         case 'V'
%                             [interpolants_V(i)] = temp;
%                             [interpolants_V_bodyframe(i)] = temp2;
%                         case 'W'
%                             [interpolants_W(i)] = temp;
%                             [interpolants_W_bodyframe(i)] = temp2;
%                     end
                    
            end
            
            
            if ( strcmp( input.kinematics_interp_method , 'trig') && isnan(temp.a0) ) || (  strcmp( input.kinematics_interp_method , 'spline') && isempty(temp)  )
                vel_field_avg_U(i) = NaN;  vel_field_avg_V(i) = NaN;  vel_field_avg_W(i) = NaN;
                vel_field_avg_U_bodyframe(i) = NaN;  vel_field_avg_V_bodyframe(i) = NaN;  vel_field_avg_W_bodyframe(i) = NaN;
            else
                switch input.kinematics_interp_method
                    case 'trig'
                        fun = @(x) trig_interp_eval(temp, x);
                        fun2 = @(x) trig_interp_eval(temp2, x);
                    case 'spline'
                        fun = @(x) fnval(temp,x);
                        fun2 = @(x) fnval(temp2,x);
                end
                
                switch names{n}
                    case 'U'
                        vel_field_avg_U(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                        vel_field_avg_U_bodyframe(i) = 1/diff([0 2*pi]) * integral(fun2, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                    case 'V'
                        vel_field_avg_V(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                        vel_field_avg_V_bodyframe(i) = 1/diff([0 2*pi]) * integral(fun2, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                    case 'W'  
                        vel_field_avg_W(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                        vel_field_avg_W_bodyframe(i) = 1/diff([0 2*pi]) * integral(fun2, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                end
            end
        end
        
    end
    
    save([outfolder,dumps{d}(1:end-8),'vel_fields_z_slice.mat']);
    return
    %%
    min_distance = 2;  do_waitbar = true;
    [OK_points] = filter_inside_or_near_boundaries(points,min_distance,Solutions,Meshes, do_waitbar);
    %%
end

return

%%


%  tail_ind = find(strcmp({Mesh.name},'Tail'));
%         body_ind = find(strcmp({Mesh.name},'Body'));
%         if ~isempty([tail_ind body_ind])
%         [body_tail] = combine_Meshes(Mesh, [body_ind tail_ind]); %%leave out transverse, wingtip since it's 2D
%         else
%             body_tail = [];
%         end
%
%         [except_body] = combine_Meshes(Mesh, setdiff(1:length(Mesh),[tail_ind body_ind]));  %this will work even if body isn't there
%         [corner_elems, flat_verts] = flatten_mesh(except_body);

c = 5;
slice_types = {'x','y','z'};
frames = {'fixed','body'};

% slice_types = {'x','y','z'};

slice_types = {'z'};
frames = {'fixed'};
%  frames = {'body'};
filter_mode = 'filtered';
do_save = false;
% do_save = true;

filter_color = false;

clear settings
settings.x.scale.fixed = 20; %40 for 2-1.5, 3.5-0    settings.x.scale.body = 4; % 20 for 2-1.5 hairs
settings.x.cutoff.fixed = 1;    settings.x.cutoff.body = 1;
settings.x.value = 30.05;
settings.x.caxis.fixed = [0 25];  settings.x.caxis.body = [0 180]; % 0 25 for 2-1.5 hairs
settings.x.scale_vec.fixed.velocity = [0 15 0]'; %10 for 2-1.5, 3.5-0    settings.x.scale_vec.body.velocity = [0 80 0]';   % 20 for 2-1.5
settings.x.scale_vec.position = [settings.x.value + 20, -40 , -40]';
settings.x.scale_vec.boxpos =  [settings.x.value + 1 , -40 , -40 - 1];
settings.x.light_position = [1 0 0];
settings.x.axes_limits = [-Inf Inf; -60 60; -45 45];

settings.y.scale.fixed = 4;     settings.y.scale.body = 1.5; %2.5 for 2-1.5
settings.y.cutoff.fixed = 1;    settings.y.cutoff.body = 1;
settings.y.value = 0;
settings.y.caxis.fixed = [0 120];  settings.y.caxis.body = [0 400];  %250 for 2-1.5
settings.y.scale_vec.fixed.velocity = [100 0 0]';     settings.y.scale_vec.body.velocity = [200 0 0]'; %150 for 2-1.5
settings.y.scale_vec.position = [-40,settings.y.value - 20 , -40]';
settings.y.scale_vec.boxpos =  [-40,  settings.y.value - 1 , -40 - 1];
settings.y.light_position = [0 -1 0];
settings.y.axes_limits = [-75 70; -Inf Inf; -60 60];

settings.z.scale.fixed = 4;      settings.z.scale.body = 1.5;
settings.z.cutoff.fixed = 1;     settings.z.cutoff.body = 1;
settings.z.value = 0;
settings.z.caxis.fixed = [0 150];  settings.z.caxis.body = [50 400]; %200 for 2-1.5
settings.z.scale_vec.fixed.velocity = [100 0 0]';     settings.z.scale_vec.body.velocity = [200 0 0]'; %150 for 2-1.5
settings.z.scale_vec.position = [-40 , -40, settings.z.value + 20]';
settings.z.scale_vec.boxpos =  [-40, -40 - 1,  settings.z.value + 1 ];
settings.z.light_position = [0 0 1];
settings.z.axes_limits = [-30 100; -100 35; -Inf Inf];

fontsize = 12;


for s = 1:length(slice_types)
    slice_type = slice_types{s};
    for f = 1:length(frames)
        frame = frames{f};
        c = c+1;
        
        
        % x slice
        scale = settings.(slice_type).scale.(frame);
        cutoff = settings.(slice_type).cutoff.(frame);
        switch frame
            case 'fixed'
                speed = sqrt(vel_field_avg_U.^2+vel_field_avg_V.^2+vel_field_avg_W.^2);
            case 'body'
                speed = sqrt(vel_field_avg_U_bodyframe.^2+vel_field_avg_V_bodyframe.^2+vel_field_avg_W_bodyframe.^2);
        end
        
        speed_cutoff = quantile(speed(:),cutoff);
        too_big = speed > speed_cutoff;
        switch frame
            case 'fixed'
                U_s = vel_field_avg_U;  V_s = vel_field_avg_V;  W_s = vel_field_avg_W;
            case 'body'
                U_s = vel_field_avg_U_bodyframe;  V_s = vel_field_avg_V_bodyframe;  W_s = vel_field_avg_W_bodyframe;
        end
        switch filter_mode
            case 'filtered'
                U_s(too_big) = NaN;  V_s(too_big) = NaN;  W_s(too_big) = NaN;
        end
        
        
        
        
        
        
        switch slice_type
            case 'x'
                x = settings.x.value;
                y = Ly(1):spacing:Ly(2);
                z = Lz(1):spacing:Lz(2);
                U_s(:) = 0;
            case 'y'
                y = settings.y.value;
                x = Lx(1):spacing:Lx(2);
                z = Lz(1):spacing:Lz(2);
                V_s(:) = 0;
            case 'z'
                z = settings.z.value;
                y = Ly(1):spacing:Ly(2);
                x = Lx(1):spacing:Lx(2);
                W_s(:) = 0;
        end
        
        Speed = sqrt( U_s.^2 + V_s.^2 + W_s.^2 );
        %         U_s = U_s*scale;  V_s = V_s*scale ;  W_s = W_s*scale;
        %
        [Xtemp,Ytemp,Ztemp] = meshgrid(x,y,z);
        %         [temp] = inside_mesh([Xtemp(:) Ytemp(:) Ztemp(:)], body_tail);
        %         is_inside = NaN(size(Xtemp));
        %         is_inside(:) = temp;  %preserve dimensions of is_inside to match X, Y, Z
        %         points = [Xtemp(~is_inside) Ytemp(~is_inside) Ztemp(~is_inside)];
        points_temp = [Xtemp(:) Ytemp(:) Ztemp(:)];
        if strcmp(filter_mode,'filtered')
            % %         leave out any internal points
            is_OK =  ismember(points_temp,OK_points,'rows');
            points_temp = points_temp(is_OK,:);
        end
        
        %         [ distances ] = point2trimesh('Faces',corner_elems,'Vertices',flat_verts,'QueryPoints',points, 'Algorithm','parallel');
        %         points = points(abs(distances) >= min_distance,:);
        
        figure(c); clf;  [s,e] = plot_mesh(Mesh); hold on;  axis tight;  lh(c) = light('position',settings.(slice_type).light_position);
        set(e,'edgealpha',0.1);
        
        switch slice_type
            case 'x'
                sl = slice(X,Y,Z,Speed,settings.x.value,[],[]);
                
            case 'y'
                sl = slice(X,Y,Z,Speed,[],settings.y.value,[]);
            case 'z'
                sl = slice(X,Y,Z,Speed,[],[],settings.z.value);
        end
        
        if strcmp(filter_mode,'filtered') && filter_color
            do_waitbar = true;
            slice_points = [sl.XData(:) sl.YData(:) sl.ZData(:)];
            [OK_slice_points] = filter_inside_or_near_boundaries(slice_points,min_distance,Solutions,Meshes, do_waitbar);
            isbad = ~ismember(slice_points, OK_slice_points,'rows');  sl.CData(isbad) = NaN;
        end
        
        
        set(sl,'facecolor','interp','edgecolor','none','facealpha',1);
        cb = colorbar;  cb.FontSize = fontsize;
        cblabel('2D speed (\mum s^{-1})','fontsize',fontsize)
        caxis(settings.(slice_type).caxis.(frame))
        hold on
        
        scale_struct.s = scale;
        scale_struct.velocity = settings.(slice_type).scale_vec.(frame).velocity;
        scale_struct.position = settings.(slice_type).scale_vec.position;
        
        % cp = coneplot(X,Y,Z,U_s,V_s,W_s,points(:,1),points(:,2),points(:,3),0,'quiver');
        %         cp = coneplot(X,Y,Z,U_s,V_s,W_s,points_temp(:,1),points_temp(:,2),points_temp(:,3),0,'quiver','nearest');
        
        [cp,svec] = coneplot(X,Y,Z,U_s,V_s,W_s,points_temp(:,1),points_temp(:,2),points_temp(:,3),scale_struct,'linear');
        alpha = 0.15;  cp.EdgeAlpha = alpha;  cp.FaceAlpha = alpha;
        
        %         cp = coneplot(X,Y,Z,U_s,V_s,W_s,points_temp(:,1),points_temp(:,2),points_temp(:,3),'nearest');
        clear text
        tx  = text(settings.(slice_type).scale_vec.boxpos(1),settings.(slice_type).scale_vec.boxpos(2),settings.(slice_type).scale_vec.boxpos(3),...
            { [num2str(norm(scale_struct.velocity)),' \mum s^{-1}'] } );
        tx.BackgroundColor = 'w'; tx.EdgeColor = 'k'; tx.HorizontalAlignment = 'center';  tx.VerticalAlignment = 'cap'; tx.LineWidth = 2; tx.Margin = 8;
        tx.FontWeight = 'bold';  tx.FontSize = fontsize;
        %         cp = QUIVER3D(X,Y,Z,U,V,W,COLOR,S,N)   may be worth trying, need
        %         to get velocities that go with points_temp
        
        
        hold on
        
        
        %         set(cp,'Color',[0 0 0]);
        % cp = coneplot(X,Y,Z,U,V,W,scale,'nointerp');
        switch slice_type
            case 'x'
                set(gca,'view',[90 0]);
                title({[frame,' frame'],'x slice through transverse sheet midpoint'},'fontsize',14)
                %                 sl = slice(X,Y,Z,Speed,settings.x.value,[],[]);
                
            case 'y'
                set(gca,'view',[0 0]);
                title({[frame,' frame'],'y slice through body center'},'fontsize',14)
                %                 sl = slice(X,Y,Z,Speed,[],settings.y.value,[]);
            case 'z'
                set(gca,'view',[0 90]);
                title({[frame,' frame'],'z slice through body center'},'fontsize',14)
                %                 sl = slice(X,Y,Z,Speed,[],[],settings.z.value);
        end
        %         set(sl,'facecolor','interp','edgecolor','none','facealpha',0.6);
        
        
        xlabel('x (\mum)');   ylabel('y (\mum)');   zlabel('z (\mum)')
        set(gca,'fontsize',fontsize);
        
        %         xlim([min(X(:)) max(X(:))]);  ylim([min(Y(:)) max(Y(:))]);  zlim([min(Z(:)) max(Z(:))]);
        xlim(settings.(slice_type).axes_limits(1,:));  ylim(settings.(slice_type).axes_limits(2,:));  zlim(settings.(slice_type).axes_limits(3,:));
        
        set(gcf,'Position',[   1          41        1600        840]);
        drawnow
        
        if do_save
            saveas(gcf,[outfolder,slice_type,' slice ',frame,' frame',' ',filter_mode,'.fig'])
            saveas(gcf,[outfolder,slice_type,' slice ',frame,' frame',' ',filter_mode,'.png'])
        end
        
    end
end
return
%% z slice  Need to run right after above cell that makes the z slice figure!
folder = 'C:\Users\rudi\Desktop\RD\exported dino fields\';

ind = find(Z(1,1,:) == settings.z.value);
Xtemp2 = flipud(X(:,:,ind))';  Ytemp2 = flipud(Y(:,:,ind))';
Utemp2 = flipud(U_s(:,:,ind))';  Vtemp2 = flipud(V_s(:,:,ind))';


inds = [];
for i = 1:numel(Xtemp2)
    if ~ ismember( round([Xtemp2(i) Ytemp2(i) settings.z.value],9) , round(points_temp,9) ,'rows' ) %point is inside body, round to avoid roundoff issues....
        Utemp2(i) = 0;  Vtemp2(i) = 0;
    end
    if ~ismember(Xtemp2(i),x) ||  ~ismember(Ytemp2(i),y)  % removes "extra" coordinates that go with slice planes
        inds(end+1) = i;
    end
end

Xtemp2(inds) = []; Ytemp2(inds) = []; Utemp2(inds) = []; Vtemp2(inds) = [];

data = [ Xtemp2(:)/1000 Ytemp2(:)/1000 Utemp2(:)*1E-6 Vtemp2(:)*1E-6 ];
data = sortrows(data, [-2 1]); % sorts first by 2nd coord descending, then within constant 2nd coord, 1st coord ascending as per DaVis requirement


fileID = fopen([folder,'Z_slice.txt'],'w');
fprintf(fileID,'#DaVis 8.0.6 2D-vector 8 %i %i "position" "mm" "position" "mm" "velocity" "m/s"\r\n',[length(unique(Ytemp2)) length(unique(Xtemp2))]);
fprintf(fileID,'%.12g\t%.12g\t%.12g\t%.12g\r\n',data');
fclose(fileID);
%% x slice  Need to run right after above cell that makes the x slice figure!
folder = 'C:\Users\rudi\Desktop\RD\exported dino fields\';

ind = find(X(1,:,1) == settings.x.value); %meshgrid, not ndgrid....
Ytemp2 = flipud(squeeze(Y(:,ind,:)))';  Ztemp2 = flipud(squeeze(Z(:,ind,:)))';
Vtemp2 = flipud(squeeze(V_s(:,ind,:)))';  Wtemp2 = flipud(squeeze(W_s(:,ind,:)))';

inds = [];
for i = 1:numel(Ytemp2)
    if ~ ismember( round([settings.x.value Ytemp2(i) Ztemp2(i) ],9) , round(points_temp,9) ,'rows' ) %point is inside body, round to avoid roundoff issues....
        Vtemp2(i) = 0;  Wtemp2(i) = 0;
    end
    if ~ismember(Ytemp2(i),y) ||  ~ismember(Ztemp2(i),z)  % removes "extra" coordinates that go with slice planes
        inds(end+1) = i;
    end
end

Ytemp2(inds) = []; Ztemp2(inds) = []; Vtemp2(inds) = []; Wtemp2(inds) = [];

data = [ Ytemp2(:)/1000 Ztemp2(:)/1000 Vtemp2(:)*1E-6 Wtemp2(:)*1E-6 ];
data = sortrows(data, [-2 1]); % sorts first by 2nd coord descending, then within constant 2nd coord, 1st coord ascending as per DaVis requirement


fileID = fopen([folder,'X_slice.txt'],'w');
fprintf(fileID,'#DaVis 8.0.6 2D-vector 8 %i %i "position" "mm" "position" "mm" "velocity" "m/s"\r\n',[length(unique(Ztemp2)) length(unique(Ytemp2))]);
fprintf(fileID,'%.12g\t%.12g\t%.12g\t%.12g\r\n',data');
fclose(fileID);
%% y slice   Need to run right after above cell that makes the y slice figure!
folder = 'C:\Users\rudi\Desktop\RD\exported dino fields\';

ind = find(Y(:,1,1) == settings.y.value); %meshgrid, not ndgrid....
Xtemp2 = flipud(squeeze(X(ind,:,:)))';  Ztemp2 = flipud(squeeze(Z(ind,:,:)))';
Utemp2 = flipud(squeeze(U_s(ind,:,:)))';  Wtemp2 = flipud(squeeze(W_s(ind,:,:)))';

inds = [];
for i = 1:numel(Xtemp2)
    if ~ ismember( round([ Xtemp2(i) settings.y.value Ztemp2(i) ],9) , round(points_temp,9) ,'rows' ) %point is inside body, round to avoid roundoff issues....
        Utemp2(i) = 0;  Wtemp2(i) = 0;
    end
    if ~ismember(Xtemp2(i),x) ||  ~ismember(Ztemp2(i),z)  % removes "extra" coordinates that go with slice planes
        inds(end+1) = i;
    end
end

Xtemp2(inds) = []; Ztemp2(inds) = []; Utemp2(inds) = []; Wtemp2(inds) = [];

data = [ Xtemp2(:)/1000 Ztemp2(:)/1000 Utemp2(:)*1E-6 Wtemp2(:)*1E-6 ];
data = sortrows(data, [-2 1]); % sorts first by 2nd coord descending, then within constant 2nd coord, 1st coord ascending as per DaVis requirement

fileID = fopen([folder,'Y_slice.txt'],'w');
fprintf(fileID,'#DaVis 8.0.6 2D-vector 8 %i %i "position" "mm" "position" "mm" "velocity" "m/s"\r\n',[length(unique(Ztemp2)) length(unique(Xtemp2))]);
fprintf(fileID,'%.12g\t%.12g\t%.12g\t%.12g\r\n',data');
fclose(fileID);