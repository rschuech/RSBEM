

if input.periodic_kinematics
    N = 2; % # initial number of equally spaced evaluation points, not including final point implicitly known to be same as first point
else
    N = 3; % # initial number of equally spaced evaluation points, including final point which is different than first point
end

phase_bounds = [0 2*pi]; %should never need to change this

fignum = input.output.interpolation.fignum;
n_pts_plot = 500;

if isfield(input,'use_interp_partial') && input.use_interp_partial && isfield(input.paths,'interp_partial_file') && exist([input.paths.dumpfolder , input.paths.interp_partial_file],'file')
    temp = load([input.paths.dumpfolder , input.paths.interp_partial_file]);
    interp_partial = temp.interp_partial;
else
    interp_partial = struct('f',{},'kinematics',{},'phase',{},'Metadata_rand_inds',{});
end


clear yprime
clear global solutions
global solutions  %global ugliness isn't necessary here, but it is for full
% timestepping simulations in which interpolation shortcut can't be used (e.g. all submeshes are deforming) and all output must be gotten from ode45, so same
% output storage method is used for both cases for consistency
% initialize phase angles here, will update at end of while loop
% solutions.phase are just the new phase angles to compute
% for.  Solutions.phase are all phase angles computed so
% far.
if input.periodic_kinematics
    solutions.phase = linspace(phase_bounds(1), phase_bounds(2),N+1);  %temporarily add last point
    solutions.phase = solutions.phase(1:end-1)'; %remove last point since it is the same as first point
else
    solutions.phase = linspace(phase_bounds(1), phase_bounds(2),N);
    solutions.phase = solutions.phase';
end


% this is going to be iter 0: the first evaluation at zero phase angle
switch input.bugtype
    case 'bacteria'
        if recycle_A %just did forced, so can recycle majority of A
           
            temp = tic;
            [ A_temp, A_motor_torque0, ~, ~, debug_info ] = matrix_assembly_mex_wrapper( Mesh,matrix_props,assembly_input );
            if input.performance.verbose
                disp(['Matrix assembly wrapper took ',num2str(toc(temp))]);
            end
            A_temp(1:matrix_props.n_col*3, 1:matrix_props.n_col*3) = A_recycle;  %this part is exactly the same as for the forced case we just did
            A_recycle = A_temp;  %stores entire freeswim A just for use during very first interpolation point, will then get deleted
        else %not recycling from the forced case, but might possibly recycle this current result for forced next

            temp = tic;
            % this is probably quite broken now due to revised velocity BC
            % handling and DL stuff
%             [A_temp, A_motor_torque0, A_force_recycle, A_torque_recycle, debug_info] = matrix_assembly_mex_wrapper(Mesh, matrix_props, assembly_input); %save A_force and A_torque for possible reuse (unlikely though)
switch input.Tail.motorBC
    case "torque"
        Mesh(1).u = zeros(Mesh(1).n_nodes,3);  Mesh(2).u = zeros(Mesh(2).n_nodes,3);
        assembly_input.Tail.motor_torque = input.Tail.motor_torque;
    case "freq"
        % need to set zero u on body, and relative u due to known freq on tail....
end
[matrix_props] = gen_matrix_props(input,Mesh, mesh_node_parameters); % this needs to be after setting BCs since we need to know about free slip nodes

% index_mapping2 = index_mapping;
% index_mapping2.local_node2global_node = cell2struct(index_mapping.local_node2global_node,'indices',1);
% index_mapping2.global_node2local_node = cell2struct(index_mapping.global_node2local_node,'indices',length(index_mapping.global_node2local_node));
t = 0;
Repulsion = calc_repulsive_forces(Mesh, Network, index_mapping, assembly_input.repulsion);

[ A_temp, A_force_recycle, A_torque_recycle, RHS, A_motor_torque0] = matrix_assembly_mex_wrapper(Mesh,Network, Repulsion, matrix_props,index_mapping,mesh_node_parameters,assembly_input, t);

            if input.performance.verbose
                disp(['Matrix assembly wrapper took ',num2str(toc(temp))]);
            end
            A_recycle = A_temp;  %and of course save A for possible reuse (unlikely though)
        end
        
        if input.performance.debug_mode
            save(input.paths.dumpfolder + input.paths.namebase.full + "_debug_info_initial" + ".mat",'debug_info');
        end
        
        %either way, now save just the parts of A that will be reused by
        %tensor rotation
        A0 = cell(length(Mesh),1);
        for ii = 1:length(Mesh)
            A0{ii} = A_temp( (Mesh(ii).indices.glob.unq_bounds.vert(1) - 1)*3+1  :  Mesh(ii).indices.glob.unq_bounds.vert(2) * 3 , ...
                (Mesh(ii).indices.glob.unq_bounds.vert(1) - 1)*3+1  :  Mesh(ii).indices.glob.unq_bounds.vert(2) * 3  );
        end
        % now A0 only contains the values it actually needs, i.e. the
        % values that will be copied and rotated between interpolation points
        % hopefully this saves RAM
        clear A_temp
        
    case {'dino', 'sheet'}
        assembly_input.skip_traction_integrals = false;
        assembly_input.bugtype = input.bugtype;
        A0 = [];  %don't bother with any fancy recycling since nothing can be recycled for dino case
        A_motor_torque0 = [];
        A_recycle = [];
end

switch input.problemtype
    case 'freeswim' %will have to expand later for dinoflagellates etc
        switch input.bugtype
            case 'bacteria'
                switch input.tail.motorBC
                    case 'freq'
                        interp_fields = {'U','Omega'};
                        labels = {'U','\Omega'};
                        isvector = [true true];
                        numplots = 6;
                    case 'torque'
                        interp_fields = {'U','Omega','omega'}; % be sure to put all scalars after all vectors
                        labels = {'U','\Omega','\omega'};
                        isvector = [true true false];
                        numplots = 7;
                end
            case {'dino', 'sheet'}
                interp_fields = {'U','Omega'};
                labels = {'U','\Omega'};
                isvector = [true true];
                numplots = 6;
        end
        nrows = any(isvector)*4 + sum(~isvector);
        ncols = sum(isvector);
end

clear Solutions
Solutions.phase = zeros(0,1);
%Solutions.f = zeros(0,matrix_props.n_col*3);
Solutions.f = cell(0,1);
Solutions.U = zeros(0,3);
Solutions.Omega = zeros(0,3);



switch input.bugtype
    case 'bacteria'
        switch assembly_input.tail.motorBC
            case 'torque'
                Solutions.omega = zeros(0,1);
            case 'freq'
                Solutions.motor_torque = zeros(0,1);
        end
end

interpolant = []; rel_diff = [];   interpolant_dir = [];  interpolant_mag = [];

switch input.bugtype
    case 'bacteria'
        if length(y0) == 6 % freq motorBC
            y0_interp = [y0; 0];  %need this so y0 matches size of y, even for freq BC when phase angle is known (just a formality for yprime() )
        else
            y0_interp = y0;
        end
    case {'dino', 'sheet'}
        y0_interp = y0;  % with fixed timestepping algorithm, y(1:3) represents a translation, which starts at zero
        Solutions.rand_inds = cell(0,1);
end

switch input.bugtype
    case 'bacteria'
        assembly_input.skip_rigid_integrals = true; %if object consists of rigidly rotating or translating submeshes (i.e. bacteria)
    case {'dino', 'sheet'}
        assembly_input.skip_rigid_integrals = false; %each dino mesh is completely different, so can't reuse any of A   :-(
end
assembly_input.skip_traction_integrals = false;  %still need to redo "cross" integrals here

t = [];  %don't need t for interpolating kinematics based on an input phase angle
n_iter = 0;

if input.output.interpolation.doplot
    x = linspace(phase_bounds(1),phase_bounds(2),n_pts_plot);
    figure(fignum)
    clf
    %     set(gcf,'position',[1163          61         751         939]);
    set(gcf,'position',[770         183        1144         789]);
    iter_text(1) = annotation('textbox',[0 0 1 0.1],'string',['adaptive iteration # ',num2str(n_iter)],'fontsize',16,'fontweight','bold','linestyle','none',...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
    
    figure(fignum+1)
    clf
    set(gcf,'position',[176    37   652   923]);
    iter_text(2) = annotation('textbox',[0 0 1 0.1],'string',['adaptive iteration # ',num2str(n_iter)],'fontsize',16,'fontweight','bold','linestyle','none',...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
    
end

while n_iter < input.accuracy.interpolation.max_iter
    n_iter = n_iter + 1;
    set(iter_text,'string',['adaptive iteration # ',num2str(n_iter)]);
    
    fields = fieldnames(Solutions);
    fields = setdiff(fields, {'phase', 'f'});  % phase angle is updated separately, f may change size for dino case
    for f = 1:length(fields)
        solutions.(fields{f}) = NaN(length(solutions.phase),size(Solutions.(fields{f}),2));
    end
    solutions.f = cell(length(solutions.phase), 1);  %each traction list will have to go into a separate cell since they may be different sizes
    solutions.n_out = 0; %reset index for last solution output
    if strcmp(input.bugtype, 'dino')
        solutions.rand_inds = cell(length(solutions.phase), 1);
    end
    
    for i = 1:length(solutions.phase)
        
        switch input.bugtype
            case 'bacteria'
                
                %                 y = [Mesh(1).refpoints(:,1); [0 0 0 solutions.phase(i)]' ];
                y = y0_interp;  y(7) = solutions.phase(i);
                Mesh_temp = Mesh;
                if ~isequal(y,y0_interp) %location and/or orientation has changed; this is not the first yprime calculation
                    
                    %first rotate only tail to handle tail phase angle relative to body
%                     switch input.tail.orientation
%                         case 'centerline'
%                             Mesh_temp(2) = rotateMesh(Mesh_temp(2),[0 0 y(7)]');  %rotate tail by y(7) around x-axis
%                         case 'pole2pole'
                              [Mesh_temp(2)] = rotateTail(Mesh_temp(2), y(7));
%                     end
                    %then move all submeshes according to translation and rotation of body +
                    %tail
                    Mesh_temp = move_Mesh(Mesh_temp,y);
                    
                end
                
            case {'dino', 'sheet'}
                if ismember(solutions.phase(i), [interp_partial.phase]) && ( ~ isfield(input,'redo_interp_phases') || ~ismember(roundn(solutions.phase(i),input.interp_phases_tol), roundn(input.redo_interp_phases,input.interp_phases_tol) ) )
                    % using interp_partial, send in existing rand_inds
                    [~,ind] = ismember(solutions.phase(i), [interp_partial.phase]);
                    [Mesh_temp, Metadata_temp, matrix_props, Metadata_rand_inds] = load_dino_mesh(solutions.phase(i), Mesh_files, input , interp_partial(ind).Metadata_rand_inds);
                else % not using interp_partial, generate new rand_inds
                    [Mesh_temp, Metadata_temp, matrix_props, Metadata_rand_inds] = load_dino_mesh(solutions.phase(i), Mesh_files, input);
                    
                end
                
                assembly_input.accuracy.include_DL = [];
                for field = string({Mesh.name})
                    if isfield(input.accuracy.include_DL,field)
                        assembly_input.accuracy.include_DL(end+1) = input.accuracy.include_DL.(field);
                    else
                        assembly_input.accuracy.include_DL(end+1) = input.accuracy.include_DL.default;
                    end
                end
                % currently the order of dino submeshes when loaded is not
                % carefully controlled so doing this here makes sure vector of eps2
                % matches order of Mesh.name
                assembly_input.accuracy.eps2 = [];
                for field = string({Mesh.name})
                    if isfield(input.accuracy.epsilon,field)
                        assembly_input.accuracy.eps2(end+1) = input.accuracy.epsilon.(field)^2;
                    else
                        assembly_input.accuracy.eps2(end+1) = input.accuracy.epsilon.default^2;
                    end
                end
                
                %                 y = [Mesh_temp(1).refpoints(:,1); [0 0 0 solutions.phase(i)]' ];
                y = y0_interp;  y(7) = solutions.phase(i);
                %                 if solutions.phase(i) == 0
                % %                     y0 = y(1:6);  %should be organized to somewhere else...
                %                      end
                if solutions.phase(i) == 0
                    Mesh = Mesh_temp; % for convenience if a mesh is needed later (e.g. for plotting) - use the initial mesh
                    % also used to get refpoint for timestepping tracking
                end
                
%                 u_rel = {};
                for n = 1:length(Mesh_temp)
                    name = Mesh_temp(n).name;
                    [~,inds] = ismember(Mesh_temp(n).indices.orig.vert,  Metadata_temp.(name).indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs in Metadata to match current order in Mesh
                    BCs.freeswim.(name) = Metadata_temp.(name).BCs(inds,:);
                    Mesh_temp(n).u = BCs.freeswim.(name); % most logical place to put relative velocities is inside Mesh ?
                    % note, for mobility problem, Mesh.u is velocity relative to
                    % body frame, whereas, for resistance problem, Mesh.u is
                    % known fixed frame velocity
                    
%                     u_rel{n} = BCs.freeswim.(name); % dynamic field references not allowed in Coder so can't use BCs struct for DL calculation
                end
                BCs.freeswim.phase_speed = Metadata_temp.geom.phase_speed;  %overall phase speed of both flagella (i.e. for a full cycle)
                input.phase_speed = Metadata_temp.geom.phase_speed;  % need to organize this at some point.  need to send value in to timestepping.m
        end
        
        
        
        % compute time derivatives for many values of tail phase angle
        %adjust tail phase angle and input to yprime():
        % note, will input tail phase angle as part of y regardless of
        % motorBC since we aren't actually timestepping yet, but imposing
        % the instantaneous configuration
        
        if input.performance.debug_mode
            debug_dump_filename = [input.paths.dumpfolder,input.paths.namebase.full,'_debug_info','_phase_',num2str(solutions.phase(i)),'.mat'];
            yprime(t,y,y0_interp,Mesh_temp,matrix_props,assembly_input,BCs,A0,A_motor_torque0, A_recycle, debug_dump_filename);
        else
            if strcmp(input.bugtype,'dino') && ismember(solutions.phase(i), [interp_partial.phase]) && ( ~ isfield(input,'redo_interp_phases') || ~ismember(roundn(solutions.phase(i),input.interp_phases_tol), roundn(input.redo_interp_phases,input.interp_phases_tol) ) )
                
                [~,ind] = ismember(solutions.phase(i), [interp_partial.phase]);
                yprime(t,y,y0_interp,Mesh_temp,matrix_props,assembly_input,BCs,A0,A_motor_torque0, A_recycle , interp_partial(ind));
                solutions.rand_inds{solutions.n_out} = Metadata_rand_inds;
                interp_partial( solutions.phase(i) == [interp_partial.phase] ).Metadata_rand_inds = Metadata_rand_inds; % no matter what, need to store rand_inds going with solution f in interp_partial
                
            else % don't have this phase in solutions_partial, need to do it from scratch
                [~,interp_partial_temp] = yprime(t,y,y0_interp,Mesh_temp,matrix_props,assembly_input,BCs,A0,A_motor_torque0, A_recycle);
                interp_partial_temp.phase = solutions.phase(i); % need to add phase here since it isn't input to yprime (though it's probably in there but7 called something else)
                solutions.rand_inds{solutions.n_out} = Metadata_rand_inds;
                interp_partial_temp.Metadata_rand_inds = Metadata_rand_inds; % no matter what, need to store rand_inds going with solution f in interp_partial
                
                ind = find(roundn(solutions.phase(i),input.interp_phases_tol) == roundn([interp_partial.phase],input.interp_phases_tol) );
                if isempty(ind) % never did this phase before
                    interp_partial(end+1) = interp_partial_temp;
                    [~,inds] = sort([interp_partial.phase]);
                    interp_partial = interp_partial(inds);  % make sure interp_partial is always sorted by phase
                else % this phase is already there but we've decided to redo and replace it
                    interp_partial(ind) = interp_partial_temp;
                    
                end
                
                save([input.paths.dumpfolder , input.paths.interp_partial_file] , 'interp_partial');
                
            end
        end
        
        if ~input.plot_discarded_interp_phases
            filter = ~ismember(roundn(solutions.phase,input.interp_phases_tol), roundn(input.discard_interp_phases,input.interp_phases_tol));
            for ff = 1:length(interp_fields)
                solutions.(interp_fields{ff})(~filter,:) = NaN;
            end
        end
        
        if recycle_A % we just did forced, now we just did very first interpolation point for freeswim and so simply copied A_recycle exactly from above
            %but now, we'll never need it again, so get rid of it to free RAM
            A_recycle = [];  %don't actually clear cause you need a placeholder to keep inputting to yprime in further loop iterations
        end %otherwise, we should keep A_recycle from above around on the off chance we do forced after freeswim
        
        if input.output.interpolation.doplot  %update plot with each new interpolation point as it's computed
            figure(fignum)
            pc = 0;
            for f = 1:length(interp_fields)
                pc = pc + 1;
                vector_mag = sqrt( sum( solutions.(interp_fields{f}).^2 , 2 ) ); %vector magnitude at each phase
                
                if isvector(f)
                    sp = f; % move along top row
                else
                    start = sum(isvector)*4 + ( f - sum(isvector) ); % go to last vector subplot, then increment by a scalar subplot counter
                    sp = [start , start + ncols - 1];
                end
                subplot(nrows,ncols,sp)
                hold on
                plot(solutions.phase(i),vector_mag(i),'go','markerfacecolor','g');
                xlim(phase_bounds)
                grid on
                set(gca,'fontsize',12);
                if isvector(f)
                    ylabel(['|',labels{f},'|'],'fontsize',16,'fontweight','bold');
                else
                    ylabel([labels{f}],'fontsize',16);
                    if f == length(interp_fields)
                        xlabel(['phase angle (rad)'],'fontsize',14);
                    end
                end
                
                for j = 1:size(solutions.(interp_fields{f}),2)
                    if ~isvector(f)
                        break % no more vector variables with direction
                    end
                    row = j + 1; col = f;  sp = col + (row - 1)*ncols;
                    subplot(nrows,ncols,sp)
                    hold on
                    plot(solutions.phase(i),solutions.(interp_fields{f})(i,j) / vector_mag(i),'go','markerfacecolor','g');
                    xlim(phase_bounds)
                    grid on
                    set(gca,'fontsize',12);
                    if all(isvector) && j == 3
                        xlabel(['phase angle (rad)'],'fontsize',14);
                    end
                    ylabel([ '$\widehat{',labels{f},'}','_{',num2str(j),'}$'],'fontsize',16,'Interpreter','latex');
                    
                end
            end
            drawnow
            %
            %             if input.output.interpolation.doprint
            %                 try
            %                     print('-dpdf','-fillpage',[input.paths.dumpfolder,input.paths.namebase.full,'_interpolation','.pdf'])
            %                     saveas(input.output.interpolation.fignum,[input.paths.dumpfolder,input.paths.namebase.full,'_interpolation','.fig']);
            %                 end
            %             end
            
        end
        
        
        
        figure(fignum+1)
        pcc = 0;
        for f = 1:length(interp_fields)
            for j = 1:size(solutions.(interp_fields{f}),2)
                pcc = pcc + 1;
                subplot(numplots,1,pcc)
                hold on
                plot(solutions.phase(i),solutions.(interp_fields{f})(i,j),'go','markerfacecolor','g');
                xlim(phase_bounds)
                grid on
                set(gca,'fontsize',12);
                
                if size(Solutions.(interp_fields{f}),2) > 1
                    ylabel([labels{f},'_',num2str(j)],'fontsize',16);
                else
                    ylabel([labels{f}],'fontsize',16);
                end
                if pcc == numplots
                    xlabel(['phase angle (rad)'],'fontsize',14);
                end
            end
        end
        drawnow
        
        
        
    end
    
    %combine old and new values into best current data set
    fields = setdiff(fieldnames(Solutions),{'f', 'rand_inds'});  %solutions may have more fields than Solutions but we only care about Solutions fields
    for f = 1:length(fields)
        temp2 = NaN(size(Solutions.(fields{f}),1)+size(solutions.(fields{f}),1),size(Solutions.(fields{f}),2)); %updated version has old + new values
        for j = 1:size(Solutions.(fields{f}),2)  %components of current variable
            if input.periodic_kinematics
                temp = [Solutions.(fields{f})(:,j)'; solutions.(fields{f})(:,j)'];  %first row is old values, 2nd row is new values
                temp2(:,j) = temp(:);  % : operator concatenates columns, so this interleaves old and newly interpolated theta values
            else
                if isempty(Solutions.(fields{f})(:,j))
                    temp2(:,j) = solutions.(fields{f})(:,j);
                else
                    temp = [Solutions.(fields{f})(1:end-1,j)'; solutions.(fields{f})(:,j)'];  %temporarily discard last value at 2*pi (from the first iter) and add back after interleaving
                    temp2(:,j) = [temp(:);  Solutions.(fields{f})(end,j)];  % : operator concatenates columns, so this interleaves old and newly interpolated theta values
                end
            end
            
        end
        Solutions.(fields{f}) = temp2;
        
    end
    
    if strcmp(input.bugtype,'dino')
        vars = {'f','rand_inds'};
    else
        vars = {'f'};
    end
    
    for v = 1:length(vars)
        var  = vars{v};
        if input.periodic_kinematics
            temp = [Solutions.(var)(:,1)'; solutions.(var)(:,1)'];  %first row is old values, 2nd row is new values
            Solutions.(var) = temp(:); %intermediate shat seems unnecessary
        else
            if isempty(Solutions.(var))
                temp = [solutions.(var)(:,1)];
                Solutions.(var) = temp(:); %intermediate shat seems unnecessary
            else
                temp = [Solutions.(var)(1:end-1,1)'; solutions.(var)(:,1)'];
                Solutions.(var) = [temp(:); Solutions.(var)(end,1)]; %intermediate shat seems unnecessary
            end
        end
    end
    % temp2(:,j) = temp(:);  % : operator concatenates columns, so this interleaves old and newly interpolated theta values
    
    
    %     if strcmp(input.bugtype,'dino')
    %         temp = [Solutions.rand_inds(:,1)'; solutions.rand_inds(:,1)'];
    %         Solutions.rand_inds = temp(:);
    %     end
    
    c = 0; current_diffs = [];
    for f = 1:length(interp_fields)
        if input.accuracy.interpolation.vector_normalization % usual vector magnitude for U and Omega.  in case of bacteria, vector magnitude of tail rotation rate omega is still omega, so no change there
            vector_mag = sqrt( sum( Solutions.(interp_fields{f}).^2 , 2 ) ); %vector magnitude at each phase
            switch input.kinematics_interp_method
                case 'trig'
                    
                    interp_mag = trig_interp_fit(Solutions.phase, vector_mag);
                case 'spline'
                    filter = ~ismember(roundn(Solutions.phase,input.interp_phases_tol), roundn(input.discard_interp_phases,input.interp_phases_tol));
                    if input.periodic_kinematics
                        interp_mag = csape([Solutions.phase(filter) ; phase_bounds(2)]' , [vector_mag(filter) ; vector_mag(1)]','periodic'); % assume that first entry is phase = 0, which is never discarded
                    else
                        interp_mag = spline([Solutions.phase(filter) ]' , [vector_mag(filter) ]');
                    end
            end
            
            interpolant_mag(n_iter).(interp_fields{f}) = interp_mag;
            
            if n_iter >= 2 %we have enough iterations to compare the last two
                switch input.kinematics_interp_method
                    case 'trig'
                        fun2 = @(x) trig_interp_eval(interpolant_mag(n_iter).(interp_fields{f}), x);   %more accurate
                        fun1 = @(x) trig_interp_eval(interpolant_mag(n_iter-1).(interp_fields{f}), x); %less accurate
                        fun_norm = @(x) trig_interp_eval(interp_mag,x);
                    case 'spline'
                        fun2 = @(x) fnval(interpolant_mag(n_iter).(interp_fields{f}), x);   %more accurate
                        fun1 = @(x) fnval(interpolant_mag(n_iter-1).(interp_fields{f}), x); %less accurate
                        fun_norm = @(x) fnval(interp_mag,x);
                end
                
                
                
                
                rel_diff.mag(n_iter).(interp_fields{f}) = RMS_fun_diff_test(fun1, fun2, fun_norm, phase_bounds);  %global phase angle should always be defined from 0 - 2*pi, and will contain entire cyclic motion
                %                 c = c +1;
                current_diffs.mag(f) = rel_diff.mag(n_iter).(interp_fields{f});
            else
                rel_diff.mag(n_iter).(interp_fields{f}) = NaN;  %placeholder for first undefined relative difference
            end
            
            
            
        end
        
        for i = 1:size(Solutions.(interp_fields{f}),2) %scalar components of this variable, e.g. U(1:3), Omega(1:3).  trig_interp only works for scalars currently
            switch input.kinematics_interp_method
                case 'trig'
                    interpolant(n_iter).(interp_fields{f})(i) = trig_interp_fit(Solutions.phase,Solutions.(interp_fields{f})(:,i));
                    interpolant_dir(n_iter).(interp_fields{f})(i) = trig_interp_fit(Solutions.phase,Solutions.(interp_fields{f})(:,i) ./ vector_mag);
                case 'spline'
                    filter = ~ismember(roundn(Solutions.phase,input.interp_phases_tol), roundn(input.discard_interp_phases,input.interp_phases_tol));
                    if input.periodic_kinematics
                        interpolant(n_iter).(interp_fields{f})(i) = csape([Solutions.phase(filter) ; phase_bounds(2)]',[Solutions.(interp_fields{f})(filter,i); Solutions.(interp_fields{f})(1,i)]','periodic');
                        interpolant_dir(n_iter).(interp_fields{f})(i) = csape([Solutions.phase(filter) ; phase_bounds(2)]',[Solutions.(interp_fields{f})(filter,i) ./ vector_mag(filter) ; Solutions.(interp_fields{f})(1,i) ./ vector_mag(1)]','periodic');
                    else
                        interpolant(n_iter).(interp_fields{f})(i) = spline([Solutions.phase(filter) ]',[Solutions.(interp_fields{f})(filter,i)]');
                        interpolant_dir(n_iter).(interp_fields{f})(i) = spline([Solutions.phase(filter) ]',[Solutions.(interp_fields{f})(filter,i) ./ vector_mag(filter) ]');
                        
                    end
                    
                    if n_iter >= 2 %we have enough iterations to compare the last two
                        switch input.kinematics_interp_method
                            case 'trig'
                                fun2 = @(x) trig_interp_eval(interpolant_dir(n_iter).(interp_fields{f})(i), x);   %more accurate
                                fun1 = @(x) trig_interp_eval(interpolant_dir(n_iter-1).(interp_fields{f})(i), x); %less accurate
                            case 'spline'
                                fun2 = @(x) fnval(interpolant_dir(n_iter).(interp_fields{f})(i), x);   %more accurate
                                fun1 = @(x) fnval(interpolant_dir(n_iter-1).(interp_fields{f})(i), x); %less accurate
                        end
                        
                        if input.accuracy.interpolation.vector_normalization
                            fun_norm = @(x) ones(size(x));
                        else
                            fun_norm = fun_2;
                        end
                        
                        rel_diff.dir(n_iter).(interp_fields{f})(i) = RMS_fun_diff_test(fun1, fun2, fun_norm, phase_bounds);  %global phase angle should always be defined from 0 - 2*pi, and will contain entire cyclic motion
                        c = c+1;
                        current_diffs.dir(c) = rel_diff.dir(n_iter).(interp_fields{f})(i);  %concatenate all rel_diffs in a vector for convenience below
                        
                        % now do original rel diffs just for giggles
                        switch input.kinematics_interp_method
                            case 'trig'
                                fun2 = @(x) trig_interp_eval(interpolant(n_iter).(interp_fields{f})(i), x);   %more accurate
                                fun1 = @(x) trig_interp_eval(interpolant(n_iter-1).(interp_fields{f})(i), x); %less accurate
                                fun_norm = @(x) trig_interp_eval( interp_mag , x);
                            case 'spline'
                                fun2 = @(x) fnval(interpolant(n_iter).(interp_fields{f})(i), x);   %more accurate
                                fun1 = @(x) fnval(interpolant(n_iter-1).(interp_fields{f})(i), x); %less accurate
                                fun_norm = @(x) fnval( interp_mag , x);
                        end
                        rel_diff.components(n_iter).(interp_fields{f})(i) = RMS_fun_diff_test(fun1, fun2, fun_norm, phase_bounds);  %global phase angle should always be defined from 0 - 2*pi, and will contain entire cyclic motion
                        %                 current_diffs.(c) = rel_diff.dir(n_iter).(interp_fields{f})(i);  %concatenate all rel_diffs in a vector for convenience below
                    else
                        rel_diff.components(n_iter).(interp_fields{f})(i) = NaN;  %placeholder for first undefined relative difference
                    end
                    
            end
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    figure(fignum)
    pc = 0;
    for f = 1:length(interp_fields)
        %                  pc = pc + 1;
        vector_mag = sqrt( sum( solutions.(interp_fields{f}).^2 , 2 ) ); %vector magnitude at each phase
        if isvector(f)
            sp = f; % move along top row
        else
            start = sum(isvector)*4 + ( f - sum(isvector) ); % go to last vector subplot, then increment by a scalar subplot counter
            sp = [start , start + ncols - 1];
        end
        subplot(nrows,ncols,sp)
        cla
        
        if n_iter >= 2 %plot 2nd to last solutions if they exist
            switch input.kinematics_interp_method
                case 'trig'
                    plot(x, trig_interp_eval(interpolant_mag(n_iter-1).(interp_fields{f}), x),'b--','linewidth',1.5);
                case 'spline'
                    plot(x, fnval(interpolant_mag(n_iter-1).(interp_fields{f}), x),'b--','linewidth',1.5);
            end
            if input.periodic_kinematics
                mag1 = sqrt(sum(  [Solutions.(interp_fields{f})(1:2:end,:); Solutions.(interp_fields{f})(1,:)].^2 , 2));
                plot([Solutions.phase(1:2:end); phase_bounds(2)],mag1 ,'bo','markerfacecolor','b');
            else
                mag1 = sqrt(sum(  [Solutions.(interp_fields{f})(1:2:end,:); ].^2 , 2));
                plot([Solutions.phase(1:2:end)],mag1 ,'bo','markerfacecolor','b');
            end
            
        end
        switch input.kinematics_interp_method
            case 'trig'
                plot(x, trig_interp_eval(interpolant_mag(n_iter).(interp_fields{f}), x),'r-','linewidth',1.5);
            case 'spline'
                plot(x, fnval(interpolant_mag(n_iter).(interp_fields{f}), x),'r-','linewidth',1.5);
        end
        if input.periodic_kinematics
            mag2 = sqrt(sum(  [Solutions.(interp_fields{f})(:,:); Solutions.(interp_fields{f})(1,:)].^2 , 2));
            plot([Solutions.phase(:); phase_bounds(2)],mag2,'ro','markerfacecolor','r');
        else
            mag2 = sqrt(sum(  [Solutions.(interp_fields{f})(:,:)].^2 , 2));
            plot([Solutions.phase(:)],mag2,'ro','markerfacecolor','r');
        end
        
        xlim(phase_bounds)
        switch input.kinematics_interp_method
            case 'trig'
                fmean = 1/diff(phase_bounds) * integral(@(x) trig_interp_eval(interpolant_mag(n_iter).(interp_fields{f}), x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12);
            case 'spline'
                fmean = 1/diff(phase_bounds) * integral(@(x) fnval(interpolant_mag(n_iter).(interp_fields{f}), x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12);
        end
        if n_iter >= 2
            title(['    relative difference = ',num2str(rel_diff.mag(n_iter).(interp_fields{f}),3),'      mean = ',num2str(fmean,3)]);
        end
        
        
        for j = 1:size(solutions.(interp_fields{f}),2)
            if ~isvector(f)
                break % no more vector variables with direction
            end
            pc = pc + 1;
            row = j + 1; col = f;  sp = col + (row - 1)*ncols;
            subplot(nrows,ncols,sp)
            %                     hold on
            cla
            if n_iter >= 2 %plot 2nd to last solutions if they exist
                switch input.kinematics_interp_method
                    case 'trig'
                        plot(x, trig_interp_eval(interpolant_dir(n_iter-1).(interp_fields{f})(j), x),'b--','linewidth',1.5);
                    case 'spline'
                        plot(x, fnval(interpolant_dir(n_iter-1).(interp_fields{f})(j), x),'b--','linewidth',1.5);
                end
                
                if input.periodic_kinematics
                    plot([Solutions.phase(1:2:end); phase_bounds(2)],[Solutions.(interp_fields{f})(1:2:end,j); Solutions.(interp_fields{f})(1,j)]./mag1 ,'bo','markerfacecolor','b');
                else
                    plot([Solutions.phase(1:2:end)],[Solutions.(interp_fields{f})(1:2:end,j)]./mag1 ,'bo','markerfacecolor','b');
                end
            end
            
            switch input.kinematics_interp_method
                case 'trig'
                    plot(x, trig_interp_eval(interpolant_dir(n_iter).(interp_fields{f})(j), x),'r-','linewidth',1.5);
                case 'spline'
                    plot(x, fnval(interpolant_dir(n_iter).(interp_fields{f})(j), x),'r-','linewidth',1.5);
            end
            if input.periodic_kinematics
                plot([Solutions.phase(:); phase_bounds(2)],[Solutions.(interp_fields{f})(:,j); Solutions.(interp_fields{f})(1,j)]./mag2,'ro','markerfacecolor','r');
            else
                plot([Solutions.phase(:)],[Solutions.(interp_fields{f})(:,j)]./mag2,'ro','markerfacecolor','r');
                
            end
            xlim(phase_bounds)
            
            grid on
            switch input.kinematics_interp_method
                case 'trig'
                    fmean = 1/diff(phase_bounds) * integral(@(x) trig_interp_eval(interpolant_dir(n_iter).(interp_fields{f})(j), x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12);
                case 'spline'
                    fmean = 1/diff(phase_bounds) * integral(@(x) fnval(interpolant_dir(n_iter).(interp_fields{f})(j), x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12);
            end
            if n_iter >= 2
                title(['    relative difference = ',num2str(rel_diff.dir(n_iter).(interp_fields{f})(j),3),'      mean = ',num2str(fmean,3)]);
                ylim([-1.5 1.5]);
            end
            
            
        end
    end
    
    
    pcc = 0;
    for f = 1:length(interp_fields)
        for i = 1:size(Solutions.(interp_fields{f}),2) %scalar components of this variable, e.g. U(1:3), Omega(1:3).  trig_interp only works for scalars currently
            
            figure(fignum+1)
            pcc = pcc+1;
            subplot(numplots,1,pcc)
            cla
            if n_iter >= 2 %plot 2nd to last solutions if they exist
                switch input.kinematics_interp_method
                    case 'trig'
                        plot(x, trig_interp_eval(interpolant(n_iter-1).(interp_fields{f})(i), x),'b--','linewidth',1.5);
                    case 'spline'
                        plot(x, fnval(interpolant(n_iter-1).(interp_fields{f})(i), x),'b--','linewidth',1.5);
                end
                if input.periodic_kinematics
                    plot([Solutions.phase(1:2:end); phase_bounds(2)],[Solutions.(interp_fields{f})(1:2:end,i); Solutions.(interp_fields{f})(1,i)] ,'bo','markerfacecolor','b');
                else
                    plot([Solutions.phase(1:2:end)],[Solutions.(interp_fields{f})(1:2:end,i)] ,'bo','markerfacecolor','b');
                end
            end
            switch input.kinematics_interp_method
                case 'trig'
                    plot(x, trig_interp_eval(interpolant(n_iter).(interp_fields{f})(i), x),'r-','linewidth',1.5);
                case 'spline'
                    plot(x,fnval(interpolant(n_iter).(interp_fields{f})(i), x),'r-','linewidth',1.5);
            end
            if input.periodic_kinematics
                plot([Solutions.phase(:); phase_bounds(2)],[Solutions.(interp_fields{f})(:,i); Solutions.(interp_fields{f})(1,i)],'ro','markerfacecolor','r');
            else
                plot([Solutions.phase(:)],[Solutions.(interp_fields{f})(:,i)],'ro','markerfacecolor','r');
                
            end
            xlim(phase_bounds)
            switch input.kinematics_interp_method
                case 'trig'
                    fmean = 1/diff(phase_bounds) * integral(@(x) trig_interp_eval(interpolant(n_iter).(interp_fields{f})(i), x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12);
                case 'spline'
                    fmean = 1/diff(phase_bounds) * integral(@(x) fnval(interpolant(n_iter).(interp_fields{f})(i), x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12);
            end
            
            if n_iter >= 2
                title(['    relative difference = ',num2str(    rel_diff.components(n_iter).(interp_fields{f})(i)  ,3  ),'      mean = ',num2str(fmean,3)]);
            end
            
            
        end
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%
    
    
    
    if input.output.interpolation.doplot && input.output.interpolation.doprint
        try
            print(fignum, '-dpdf','-fillpage',[input.paths.dumpfolder,input.paths.namebase.full,'_interpolation','.pdf'])
            saveas(fignum,[input.paths.dumpfolder,input.paths.namebase.full,'_interpolation','.fig']);
            
            print(fignum+1,'-dpdf',[input.paths.dumpfolder,input.paths.namebase.full,'_interpolation_components','.pdf'])
            saveas(fignum+1,[input.paths.dumpfolder,input.paths.namebase.full,'_interpolation_components','.fig']);
            
        end
    end
    
    
    if strcmp(input.bugtype,'dino') %save a dump each iteration so intermediate results can be opened elsewhere as run continues
        temp = who;
        temp = setdiff(temp,{'A_recycle','A_force_recycle','A_torque_recycle', 'Inputs', 'A0', 'A', 'debug_info'}); %don't want to save A_recycle in the .mat file
        allvars = whos;
        tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
        tosave = {allvars(tosave).name};
        temp = intersect(temp,tosave);
        save([input.paths.dumpfolder,input.paths.namebase.full,'_dump','.mat'],temp{:});
    end
    
    
    if ~isempty(current_diffs) && max(current_diffs.mag) < input.accuracy.interpolation.max_rel_diff.mag && max(current_diffs.dir) < input.accuracy.interpolation.max_rel_diff.dir  %we've allegedly achieved convergence
        %         if strcmp(input.bugtype,'bacteria') %for bacteria, only bother
        %         saving final output since the runs go relatively quickly
        %         (actually dump gets saved again anyway in main.m)
        %         save([input.paths.dumpfolder,input.paths.namebase.full,'_dump','.mat']);
        %         end
        break
    else %refine grid of phase angles
        % solutions.phase are just the new phase angles to compute
        % for.  Solutions.phase are all phase angles computed so
        % far.
        if input.periodic_kinematics
            temp = [Solutions.phase; phase_bounds(2)]; %values at phase_bounds(2) not stored since they are same as at phase_bounds(1), but need it for next line
        else
            temp = Solutions.phase;
        end
        solutions.phase = ( temp(1:end-1) + temp(2:end) ) / 2; %new theta values halfway in between old values (grid spacing is refined by factor of 2 each iteration)
        
    end
    
end %adaptive while loop
