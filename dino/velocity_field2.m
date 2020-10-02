
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
folder = 'C:\Users\rudi\Desktop\RD\dino_dumps_hairs_2_1.5_phases\';
dumps = dir([folder,'*dump.mat']);
dumps = {dumps.name};

for d = 4 %1:length(dumps)
    [d length(dumps)]
    dump = dumps{d}
    
    
    load([folder,dump],'Mesh_files','Solutions','input');
    
    
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
    
    spacing = 2.6;  %put vectors every this many microns
    
    Lx = [-75 70];
    Ly = [-65 65];
    Lz = [-65 65];
    
    
    
    
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
        Meshes{phase_ind} = Mesh;
        
        
        
        %     [body_tail] = combine_Meshes(Mesh, [1 2]); %%leave out transverse, wingtip since it's 2D
        
        % [body_tail] = combine_Meshes(Mesh, [1 ]);
        
        % points inside actually change over phase due to points inside tail, so
        % just calc velocities everywhere, why not
        %     [temp] = inside_mesh([Xlist(:) Ylist(:) Zlist(:)], body_tail);
        %     is_inside = NaN(size(Xlist));
        %     is_inside(:) = temp;  %preserve dimensions of is_inside to match X, Y, Z
        %     points = [Xlist(~is_inside) Ylist(~is_inside) Zlist(~is_inside)];
        
        
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
        r_cross_omega = cross( points - repmat(Mesh(1).refpoints(:,1)',size(points,1),1) ,  repmat(Solutions.Omega(phase_ind,:),size(points,1),1)) ;
        temp = field_vel - repmat(Solutions.U(phase_ind,:),size(field_vel,1),1) - r_cross_omega;
        
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
    
    interpolants_U = struct('x0',cell(size(vel_fields(1).X)));  interpolants_U_bodyframe = interpolants_U;
    interpolants_V = struct('x0',cell(size(vel_fields(1).X)));  interpolants_V_bodyframe = interpolants_V;
    interpolants_W = struct('x0',cell(size(vel_fields(1).X)));  interpolants_W_bodyframe = interpolants_W;
    names = {'U','V','W'};
    vel_field_avg_U = NaN(size(U));  vel_field_avg_V = NaN(size(V));  vel_field_avg_W = NaN(size(W));
    vel_field_avg_U_bodyframe = NaN(size(U));  vel_field_avg_V_bodyframe = NaN(size(V));  vel_field_avg_W_bodyframe = NaN(size(W));
    
    parfor i = 1:numel(vel_fields(1).X)
        for n = 1:length(names)
            j = [];
            y = []; yb = [];
            for j = 1:length(Solutions.phase)
                y(j) = vel_fields(j).(names{n})(i);
                yb(j) = vel_fields(j).bodyframe.(names{n})(i);
            end
            %          if any(isnan(y))
            % %             vel_field_avg_U(i) = NaN;  vel_field_avg_V(i) = NaN;  vel_field_avg_W(i) = NaN;
            %           ffffff
            % %             continue
            %         end
            temp = trig_interp_fit(Solutions.phase,y);
            temp2 = trig_interp_fit(Solutions.phase,yb);
            fields = fieldnames(temp);
            for f = 1:length(fields)
                switch names{n}
                    case 'U'
                        [interpolants_U(i).(fields{f})] = temp.(fields{f});
                        [interpolants_U_bodyframe(i).(fields{f})] = temp2.(fields{f});
                        
                    case 'V'
                        [interpolants_V(i).(fields{f})] = temp.(fields{f});
                        [interpolants_V_bodyframe(i).(fields{f})] = temp2.(fields{f});
                    case 'W'
                        [interpolants_W(i).(fields{f})] = temp.(fields{f});
                        [interpolants_W_bodyframe(i).(fields{f})] = temp2.(fields{f});
                end
            end
            
            if isnan(temp.a0)
                vel_field_avg_U(i) = NaN;  vel_field_avg_V(i) = NaN;  vel_field_avg_W(i) = NaN;
                vel_field_avg_U_bodyframe(i) = NaN;  vel_field_avg_V_bodyframe(i) = NaN;  vel_field_avg_W_bodyframe(i) = NaN;
            else
                
                switch names{n}
                    case 'U'
                        
                        fun = @(x) trig_interp_eval(temp, x);
                        vel_field_avg_U(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                        fun = @(x) trig_interp_eval(temp2, x);
                        vel_field_avg_U_bodyframe(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                        
                    case 'V'
                        fun = @(x) trig_interp_eval(temp, x);
                        vel_field_avg_V(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                        fun = @(x) trig_interp_eval(temp2, x);
                        vel_field_avg_V_bodyframe(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                    case 'W'
                        fun = @(x) trig_interp_eval(temp, x);
                        vel_field_avg_W(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                        fun = @(x) trig_interp_eval(temp2, x);
                        vel_field_avg_W_bodyframe(i) = 1/diff([0 2*pi]) * integral(fun, 0,2*pi,'reltol',1E-9,'abstol',1E-12);
                end
            end
        end
        
    end
    
    save([folder,dumps{d}(1:end-8),'vel_fields.mat']);
    
    
    %% Filter out points inside body or tail and points within distance tolerance of moving boundaries
    
    min_distance = 2.5; % 3
    
    OK_points = points;
    fh = waitbar(0,'points filtering');
    for phase_ind = 1:length(Solutions.phase)
        waitbar(phase_ind / length(Solutions.phase) );
        Mesh = Meshes{phase_ind};
        
        tail_ind = find(strcmp({Mesh.name},'Tail'));
        body_ind = find(strcmp({Mesh.name},'Body'));
        if ~isempty([tail_ind body_ind])
            [body_tail] = combine_Meshes(Mesh, [body_ind tail_ind]); %%leave out transverse, wingtip since it's 2D
            [is_inside] = inside_mesh(OK_points, body_tail);
        else
            is_inside = false(size(OK_points,1),1);
        end
        
        [except_body] = combine_Meshes(Mesh, setdiff(1:length(Mesh),[tail_ind body_ind]));  %this will work even if body isn't there
        
        [corner_elems, flat_verts] = flatten_mesh(except_body);
        [ distances ] = point2trimesh('Faces',corner_elems,'Vertices',flat_verts,'QueryPoints',OK_points, 'Algorithm','parallel');
        too_close = abs(distances) <= min_distance;
        bad_points = (is_inside | too_close);
        
        OK_points(bad_points,:) = [];
        
        
    end
    close(fh);
    
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


slice_types = {'x','y','z'};
frames = {'fixed','body'};

slice_types = {'x','y','z'};

slice_types = {'x'};
frames = {'fixed'};
% frames = {'body'};
filter_mode = 'filtered';
do_save = false;
% do_save = true;



clear settings
settings.x.scale.fixed = 1;    % settings.x.scale.fixed = 0.1;
settings.x.scale.body = 0.05  *  3;
settings.x.cutoff.fixed = 1;  %settings.x.cutoff.fixed = 0.93;  % settings.x.cutoff.fixed = 1;
settings.x.cutoff.body = 1;
settings.x.value = 30.05;

settings.y.scale.fixed = 0.1;
settings.y.scale.body = 0.03;  settings.y.scale.body = 0.015;
settings.y.cutoff.fixed = 0.9; % settings.y.cutoff.fixed = 1;
settings.y.cutoff.body = 0.999;  settings.y.cutoff.body = 1;
settings.y.value = 0;

settings.z.scale.fixed = 0.1;
settings.z.scale.body = 0.03;
settings.z.cutoff.fixed = 0.9; % settings.z.cutoff.fixed = 1;
settings.z.cutoff.body = 0.999;  settings.z.cutoff.body = 1;
settings.z.value = 0;

c = 43;
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
        
        
        U_s(:) = 0;
        
        Speed = sqrt( U_s.^2 + V_s.^2 + W_s.^2 );
        U_s = U_s*scale;  V_s = V_s*scale ;  W_s = W_s*scale;
        
        switch slice_type
            case 'x'
                x = settings.x.value;
                y = Ly(1):spacing:Ly(2);
                z = Lz(1):spacing:Lz(2);
            case 'y'
                y = settings.y.value;
                x = Lx(1):spacing:Lx(2);
                z = Lz(1):spacing:Lz(2);
            case 'z'
                z = settings.z.value;
                y = Ly(1):spacing:Ly(2);
                x = Lx(1):spacing:Lx(2);
        end
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
        
        
        
        
        figure(c); clf;  [s,e] = plot_mesh(Mesh); hold on;  axis tight;  light;
        set(e,'edgealpha',0.1);
        
        switch slice_type
            case 'x'
                sl = slice(X,Y,Z,Speed,settings.x.value,[],[]);
%                 u_sl = slice(X,Y,Z,U_s,settings.x.value,[],[]);
%                 v_sl = slice(X,Y,Z,V_s,settings.x.value,[],[]);
%                 w_sl = slice(X,Y,Z,W_s,settings.x.value,[],[]);
                
            case 'y'
                sl = slice(X,Y,Z,Speed,[],settings.y.value,[]);
            case 'z'
                sl = slice(X,Y,Z,Speed,[],[],settings.z.value);
        end
        
        set(sl,'facecolor','interp','edgecolor','none','facealpha',1);
        cb = colorbar;
        caxis([5 30]);
        hold on
        
        scale_struct.s = 40;
        scale_struct.velocity = [0 10 0]';  scale_struct.position = [settings.x.value + 20, -40 , -40]';
        
        % cp = coneplot(X,Y,Z,U_s,V_s,W_s,points(:,1),points(:,2),points(:,3),0,'quiver');
        %         cp = coneplot(X,Y,Z,U_s,V_s,W_s,points_temp(:,1),points_temp(:,2),points_temp(:,3),0,'quiver','nearest');
        
        [cp,svec] = coneplot(X,Y,Z,U_s,V_s,W_s,points_temp(:,1),points_temp(:,2),points_temp(:,3),scale_struct,'linear');
        alpha = 0.15;  cp.EdgeAlpha = alpha;  cp.FaceAlpha = alpha;
        
        %         cp = coneplot(X,Y,Z,U_s,V_s,W_s,points_temp(:,1),points_temp(:,2),points_temp(:,3),'nearest');
        
        tx  = text(settings.x.value + 1,-40,-40 - 1,{ [num2str(norm(scale_struct.velocity)),' \mum s^{-1}'] } );   
        tx.BackgroundColor = 'w'; tx.EdgeColor = 'k'; tx.HorizontalAlignment = 'center';  tx.VerticalAlignment = 'cap'; tx.LineWidth = 2; tx.Margin = 8;

        %         cp = QUIVER3D(X,Y,Z,U,V,W,COLOR,S,N)   may be worth trying, need
%         to get velocities that go with points_temp
        
        
        hold on
        
        
%         set(cp,'Color',[0 0 0]);
        % cp = coneplot(X,Y,Z,U,V,W,scale,'nointerp');
        switch slice_type
            case 'x'
                set(gca,'view',[90 0]);
                title({[frame,' frame'],'x slice through transverse groove (scale bar = 100 \mum s^{-1})'},'fontsize',14)
                %                 sl = slice(X,Y,Z,Speed,settings.x.value,[],[]);
                
            case 'y'
                set(gca,'view',[0 0]);
                title({[frame,' frame'],'y slice through body center (scale bar = 100 \mum s^{-1})'},'fontsize',14)
                %                 sl = slice(X,Y,Z,Speed,[],settings.y.value,[]);
            case 'z'
                set(gca,'view',[0 90]);
                title({[frame,' frame'],'z slice through body center (scale bar = 100 \mum s^{-1})'},'fontsize',14)
                %                 sl = slice(X,Y,Z,Speed,[],[],settings.z.value);
        end
        %         set(sl,'facecolor','interp','edgecolor','none','facealpha',0.6);
        
        pauseafra
        %get the data from regular quiver
        Ua = cp.UData;
        Va = cp.VData;
        Wa = cp.WData;
        Xa = cp.XData;
        Ya = cp.YData;
        Za = cp.ZData;
        delete(cp);
        
        arrowscale = 3  /  4;
        start = [Xa Ya Za];
        stop = start + arrowscale*[Ua Va Wa];
        axis(axis)
        
        arrow_length = 7;
        
        [h] = arrow(start(~isnan(Va),:), stop(~isnan(Va),:),'length',arrow_length,'width',1,'tipangle',40,'baseangle',35);  %length was 6
        
        xlim([min(X(:)) max(X(:))]);  ylim([min(Y(:)) max(Y(:))]);  zlim([min(Z(:)) max(Z(:))]);
        
        
        switch slice_type
            case 'x'
                coord = [x -40 -40];  vel = [0 100*scale 0];  n = 3;
            case 'y'
                coord = [-50 y -50];  vel = [ 100*scale 0 0];  n = 3;
            case 'z'
                coord = [-50 -50 z];  vel = [ 100*scale 0 0];  n = 3;
        end
        xr = [coord(1)-1 coord(1) coord(1)+1]; yr = [coord(2)-1 coord(2) coord(2)+1]; zr = [coord(3)-1 coord(3) coord(3)+1];
        [Xr,Yr,Zr] = meshgrid(xr,yr,zr);
        
        % scale arrow
        cpr = coneplot(Xr,Yr,Zr,repmat(vel(1),n,n,n),repmat(vel(2),n,n,n),repmat(vel(3),n,n,n),coord(1),coord(2),coord(3),0,'quiver');
        Uaa = cpr.UData;
        Vaa = cpr.VData;
        Waa = cpr.WData;
        Xaa = cpr.XData;
        Yaa = cpr.YData;
        Zaa = cpr.ZData;
        delete(cpr);
        
        
        start = [Xaa Yaa Zaa];
        stop = start + arrowscale*[Uaa Vaa Waa];
        axis(axis)
        
        
        [hr] = arrow(start(~isnan(Vaa),:), stop(~isnan(Vaa),:),'length',arrow_length,'width',1.5);
        % set(hr,'LineWidth',0.5);
        xlabel('x (\mum)');   ylabel('y (\mum)');   zlabel('z (\mum)')
        set(gca,'fontsize',12);
        
        xlim([min(X(:)) max(X(:))]);  ylim([min(Y(:)) max(Y(:))]);  zlim([min(Z(:)) max(Z(:))]);
        set(gcf,'Position',[   1          41        1800        1040]);
        drawnow
        if do_save
            saveas(gcf,['C:\Users\rudi\Desktop\RD\dino velocity fields\',slice_type,' slice ',frame,' frame',' ',filter_mode,'.fig'])
            saveas(gcf,['C:\Users\rudi\Desktop\RD\dino velocity fields\',slice_type,' slice ',frame,' frame',' ',filter_mode,'.png'])
        end
    end
end
return

%%
% middle of groove
scale = 0.15   /20;
speed = sqrt(vel_field_avg_U.^2+vel_field_avg_V.^2+vel_field_avg_W.^2);
cutoff = 0.995;  cutoff = 1;
speed_cutoff = quantile(speed(:),cutoff);
too_big = speed > speed_cutoff;
U_s = vel_field_avg_U;  V_s = vel_field_avg_V;  W_s = vel_field_avg_W;
U_s(too_big) = NaN;  V_s(too_big) = NaN;  W_s(too_big) = NaN;
U_s = U_s*scale;  V_s = V_s*scale ;  W_s = W_s*scale;
% U_s = U*scale;  V_s = V*scale ;  W_s = W*scale;


y = 0; %top of groove % 12.5
x = Lx(1):spacing:Lx(2);
z = Lz(1):spacing:Lz(2);
[Xtemp,Ytemp,Ztemp] = meshgrid(x,y,z);
% [temp] = inside_mesh([Xtemp(:) Ytemp(:) Ztemp(:)], body_tail);
% is_inside = NaN(size(Xtemp));
% is_inside(:) = temp;  %preserve dimensions of is_inside to match X, Y, Z
% points = [Xtemp(~is_inside) Ytemp(~is_inside) Ztemp(~is_inside)];
points = [Xtemp(:) Ytemp(:) Ztemp(:)];

figure(995);  [s,e] = plot_mesh(Mesh); hold on;  axis tight;  light;
cp = coneplot(X,Y,Z,U_s,V_s,W_s,points(:,1),points(:,2),points(:,3),0,'quiver');


% cp = coneplot(X,Y,Z,U,V,W,scale,'nointerp');
title('middle of groove (scale bar = 100 \mum s^{-1})')
set(gca,'view',[0 0]);
hold on



%get the data from regular quiver
Ua = cp.UData;
Va = cp.VData;
Wa = cp.WData;
Xa = cp.XData;
Ya = cp.YData;
Za = cp.ZData;
delete(cp);

arrowscale = 3/4;
start = [Xa Ya Za];
stop = start + arrowscale*[Ua Va Wa];
axis(axis)


[h] = arrow(start(~isnan(Va),:), stop(~isnan(Va),:),'length',8);

xlim([min(X(:)) max(X(:))]);  ylim([min(Y(:)) max(Y(:))]);  zlim([min(Z(:)) max(Z(:))]);


coord = [x -40 -40];  vel = [0 100*scale 0];  n = 3;
xr = [coord(1)-1 coord(1) coord(1)+1]; yr = [coord(2)-1 coord(2) coord(2)+1]; zr = [coord(3)-1 coord(3) coord(3)+1];
[Xr,Yr,Zr] = meshgrid(xr,yr,zr);
% scale arrow
cpr = coneplot(Xr,Yr,Zr,repmat(vel(1),n,n,n),repmat(vel(2),n,n,n),repmat(vel(3),n,n,n),coord(1),coord(2),coord(3),0,'quiver');

Uaa = cpr.UData;
Vaa = cpr.VData;
Waa = cpr.WData;
Xaa = cpr.XData;
Yaa = cpr.YData;
Zaa = cpr.ZData;
delete(cpr);


start = [Xaa Yaa Zaa];
stop = start + arrowscale*[Uaa Vaa Waa];
axis(axis)


[hr] = arrow(start(~isnan(Vaa),:), stop(~isnan(Vaa),:),'length',8);
% set(hr,'LineWidth',0.5);


xlim([min(X(:)) max(X(:))]);  ylim([min(Y(:)) max(Y(:))]);  zlim([min(Z(:)) max(Z(:))]);

%%
% horizontal thru body
scale = 0.25   /10;
speed = sqrt(vel_field_avg_U.^2+vel_field_avg_V.^2+vel_field_avg_W.^2);
cutoff = 0.975; % cutoff = 1;
speed_cutoff = quantile(speed(:),cutoff);
too_big = speed > speed_cutoff;
U_s = vel_field_avg_U;  V_s = vel_field_avg_V;  W_s = vel_field_avg_W;
U_s(too_big) = NaN;  V_s(too_big) = NaN;  W_s(too_big) = NaN;
U_s = U_s*scale;  V_s = V_s*scale ;  W_s = W_s*scale;
% U_s = U*scale;  V_s = V*scale ;  W_s = W*scale;


z = 0; %top of groove % 12.5
y = Ly(1):spacing:Ly(2);
x = Lx(1):spacing:Lx(2);
[Xtemp,Ytemp,Ztemp] = meshgrid(x,y,z);
% [temp] = inside_mesh([Xtemp(:) Ytemp(:) Ztemp(:)], body_tail);
% is_inside = NaN(size(Xtemp));
% is_inside(:) = temp;  %preserve dimensions of is_inside to match X, Y, Z
% points = [Xtemp(~is_inside) Ytemp(~is_inside) Ztemp(~is_inside)];
points = [Xtemp(:) Ytemp(:) Ztemp(:)];

figure(1000);  [s,e] = plot_mesh(Mesh); hold on;  axis tight;  light;
cp = coneplot(X,Y,Z,U_s,V_s,W_s,points(:,1),points(:,2),points(:,3),0,'quiver');

set(cp,'Color',[0 0 0]);
title('middle of body horizontal (scale bar = 100 \mum s^{-1})')

set(gca,'view',[0 90]);
hold on



%get the data from regular quiver
Ua = cp.UData;
Va = cp.VData;
Wa = cp.WData;
Xa = cp.XData;
Ya = cp.YData;
Za = cp.ZData;
delete(cp);

arrowscale = 3/4;
start = [Xa Ya Za];
stop = start + arrowscale*[Ua Va Wa];
axis(axis)


[h] = arrow(start(~isnan(Va),:), stop(~isnan(Va),:),'length',7);

xlim([min(X(:)) max(X(:))]);  ylim([min(Y(:)) max(Y(:))]);  zlim([min(Z(:)) max(Z(:))]);


coord = [59 -60 z];  vel = [0 100*scale 0];  n = 3;
xr = [coord(1)-1 coord(1) coord(1)+1]; yr = [coord(2)-1 coord(2) coord(2)+1]; zr = [coord(3)-1 coord(3) coord(3)+1];
[Xr,Yr,Zr] = meshgrid(xr,yr,zr);
% scale arrow
cpr = coneplot(Xr,Yr,Zr,repmat(vel(1),n,n,n),repmat(vel(2),n,n,n),repmat(vel(3),n,n,n),coord(1),coord(2),coord(3),0,'quiver');

Uaa = cpr.UData;
Vaa = cpr.VData;
Waa = cpr.WData;
Xaa = cpr.XData;
Yaa = cpr.YData;
Zaa = cpr.ZData;
delete(cpr);


start = [Xaa Yaa Zaa];
stop = start + arrowscale*[Uaa Vaa Waa];
axis(axis)


[hr] = arrow(start(~isnan(Vaa),:), stop(~isnan(Vaa),:),'length',7);
% set(hr,'LineWidth',0.5);


xlim([min(X(:)) max(X(:))]);  ylim([min(Y(:)) max(Y(:))]);  zlim([min(Z(:)) max(Z(:))]);