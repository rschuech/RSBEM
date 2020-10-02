                                                                                                          %   long-range  short-range
% dump_folder = 'C:\Users\rudi\Desktop\RD\base_case_2_1.5\';
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_coplanar_3_normal_1.5\'; %   0.006572    0.0066685
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_tail_transverse_coplanar_3\'; %   0.0080387    0.0080848
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_centered_sheet_tail\';                   %   0.0021083   0.002114
dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_coplanar_2_normal_1.5\'; %   0.0025337   0.0050506  *
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse\';                            %   0.0036085   0.0036112
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_coplanar_2_normal_1.5\';      %   0.0024933   0.0049395  *
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_coplanar_2.5\';          %   0.006954    0.0078203  -
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_normal_1.5\';            %   0.0010599   0.0012842
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail\';                       %   0.0028131   0.0034168  -
%  dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_centered_tail\';                        %   0.0013823   0.0013865
%  dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_centered_fat_tail\';                    %   0.0013733   0.0013729
%  dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_centered_fat_sheet_tail\';              %   0.0023831   0.0024163
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_tail\';                                  %   NaN         0.00067229
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\transverse\';                                 %   0.02139     0.021462
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_interactions_off\';           %   0.0082048   0.0082243
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_coplanar_2_normal_1.5_interactions_off\'; % 0.0014964   0.001518
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\tail\';                                       %   0.0069034   0.007208
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_coplanar_2\';                            %   0.0060864   0.0065921
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_coplanar_2\';                 %   0.0068022   0.0072686
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_coplanar_2\';            %   0.005103    0.007124  -
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\normal_1.5\';                                 %   0.005103    0.007124  -
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\coplanar_2_normal_1.5\';                      %   0.0051048   0.0071366

% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\coplanar_3\';                                   %   0.0069404   0.0070865
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\coplanar_2.5\';                                 %   0.0069285   0.007381
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\coplanar_2.25\';                                 %   0.0063063   0.0074027  *
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\coplanar_2.125\';                                 %   0.005107  0.0073776 *
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\coplanar_2.0625\';                                 %   0.0039024   0.0073869  *
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\coplanar_2\';                                 %   0.00182   0.007374  *
% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\coplanar_1.75\';                                 %   0.0050805   0.0072946  *

% dump_folder = 'C:\Users\rudi\Desktop\RD\pape main results\body_transverse_tail_coplanar_3_normal_1.5_again\';



dumps = dir([dump_folder,'*dump.mat']);
dumps = {dumps.name};

timestepping_dumps = dir([dump_folder,'*_timestepping.mat']);
timestepping_dumps = {timestepping_dumps.name};

if length(dumps) > 1 || length(timestepping_dumps) > 1
    stopafra
end

for d = 1 %1:length(dumps)
    [d length(dumps)];
    dump = dumps{d}
    
    dont_load = {'dump_folder','dumps','d','timestepping_dumps'};
    varlist =strjoin(dont_load','$|'); %Join into string, separating vars by '|'
    load([dump_folder,dump],'-regexp', ['^(?!' varlist ')\w']);  % orig list was  'Mesh_files','Solutions','input'
    
    load([dump_folder,timestepping_dumps{d}],'fits');

    %%
    Meshes = [];
    powers = [];
    % average over a phase cycle
    %     fh = waitbar(0,'phase cycle velocity field avg');
    for phase_ind = 1:length(Solutions.phase)
        %         waitbar(phase_ind / length(Solutions.phase) );
        disp(['On phase ind ',num2str(phase_ind),' of ',num2str(length(Solutions.phase))]);
        % phase_ind = 1;
        phase = Solutions.phase(phase_ind);
        
        [Mesh, Metadata, matrix_props, Metadata_rand_inds] = load_dino_mesh(phase, Mesh_files, input, Solutions.rand_inds{phase_ind});
        
        Meshes{phase_ind} = Mesh;
        temp = {Mesh.name};
        
        clear BCs verts levers angular_vels tot_vels BC_vels
        
        switch input.bugtype
            case 'dino'
                for n = 1:length(Mesh)
                    name = Mesh(n).name;
                    
                    [~,inds] = ismember(Mesh(n).indices.orig.vert,  Metadata.(name).indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs in Metadata to match current order in Mesh
                    BCs.freeswim.(name) = Metadata.(name).BCs(inds,:);
                    
                    verts.(name) = Mesh(n).verts(inds,:);
                    levers.(name) = Mesh(n).verts(inds,:) - repmat(Mesh(n).refpoints,1,Mesh(n).n_vert)';
                    angular_vels.(name) = cross(  repmat(Solutions.Omega(phase_ind,:),Mesh(n).n_vert,1) , levers.(name) , 2);
                    tot_vels.(name) = repmat(Solutions.U(phase_ind,:),Mesh(n).n_vert,1) + angular_vels.(name);
                end
                
            case 'bacteria'
                % body frame velocities
                BC_vels = zeros(matrix_props.n_col*3,1);
                r = Mesh(2).verts - repmat(Mesh(2).refpoints(:,1)',Mesh(2).n_vert,1);  %use tail refpoint since it's always on motor axis
                BC_vels((Mesh(2).indices.glob.unq_bounds.vert(1)-1)*3+1  :3:  Mesh(2).indices.glob.unq_bounds.vert(2)*3,1) = Solutions.omega(phase_ind)*(Mesh(1).orientation(2)*r(:,3) - Mesh(1).orientation(3)*r(:,2));  %x equations
                BC_vels((Mesh(2).indices.glob.unq_bounds.vert(1)-1)*3+2  :3:  Mesh(2).indices.glob.unq_bounds.vert(2)*3,1) = Solutions.omega(phase_ind)*(Mesh(1).orientation(3)*r(:,1) - Mesh(1).orientation(1)*r(:,3));  %y equations
                BC_vels((Mesh(2).indices.glob.unq_bounds.vert(1)-1)*3+3  :3:  Mesh(2).indices.glob.unq_bounds.vert(2)*3,1) = Solutions.omega(phase_ind)*(Mesh(1).orientation(1)*r(:,2) - Mesh(1).orientation(2)*r(:,1));  %z equations
                
                for n = 1:length(Mesh)
                    levers{n} = Mesh(n).verts - repmat(Mesh(2).refpoints,1,Mesh(n).n_vert)';
                    angular_vels{n} = cross(  repmat(Solutions.Omega(phase_ind,:),Mesh(n).n_vert,1) , levers{n} , 2);
                    tot_vels{n} = repmat(Solutions.U(phase_ind,:),Mesh(n).n_vert,1) + angular_vels{n};
                end
        end
        
        
        
        tot_vels_vec = zeros(matrix_props.n_col*3,1);
        temp = [] ;  %x, y, z velocity at each unq vert
        for i = 1:length(Mesh)
            inds = find( Mesh(i).indices.glob.vert >= Mesh(i).indices.glob.unq_bounds.vert(1) &   Mesh(i).indices.glob.vert <= Mesh(i).indices.glob.unq_bounds.vert(2)); %local inds of "global" verts belonging to this mesh
            switch input.bugtype
                case 'dino'
                    temp = [temp; tot_vels.(Mesh(i).name)(inds,:)];
                case 'bacteria'
                    temp = [temp; tot_vels{i}(inds,:)];
            end
        end
        
        tot_vels_vec(1:3:matrix_props.n_col*3) = temp(:,1);  %slot in x-velocities
        tot_vels_vec(2:3:matrix_props.n_col*3) = temp(:,2);  %slot in y-velocities
        tot_vels_vec(3:3:matrix_props.n_col*3) = temp(:,3);  %slot in z-velocities
        
        switch input.bugtype
            case 'dino'
                [RHS] = assemble_RHS(input,Mesh, matrix_props, BCs);
                BC_vels = RHS(1:end-6); % last 6 values of RHS are for kinematic constraints, not velocity BCs
        end
        
        fmat = (reshape( Solutions.f{phase_ind} ,3,[]));  %row is x y z traction, col is global vert
        velmat = (reshape( BC_vels + tot_vels_vec ,3,[])); % add BCs in body frame to velocities due to moving body to obtain fixed frame vel at each vert
        
        f_dot_v = sum( (fmat .* velmat) , 1)';
        
        clear traction velocity fdotv
        for sf = 1:length(Mesh)  %body transverse tail wingtip
            glob_inds = Mesh(sf).indices.glob.vert;
            traction{sf} = fmat(:,glob_inds);
            velocity{sf} = velmat(:,glob_inds);
            fdotv{sf} = abs( f_dot_v(glob_inds) );
        end
%         
%%   
pause
figure(334); clf
ind = 6;  var = traction{ind};  
mag = sqrt(sum(var.^2,1));  cutoff = quantile(mag, 1);
var(:,mag > cutoff) = NaN;
[s,e] = plot_mesh(Mesh(ind),3,fdotv(ind),[0 0.7E-3]);  set(s,'facealpha',0.7); %set(e,'edgealpha',0.1);
axis tight;  view([-30 -40]);
hold on
qt = quiver3(Mesh(ind).verts(:,1),Mesh(ind).verts(:,2),Mesh(ind).verts(:,3),var(1,:)',var(2,:)',var(3,:)',1,'r','linewidth',2);
% qt = coneplot(Mesh(ind).verts(:,1),Mesh(ind).verts(:,2),Mesh(ind).verts(:,3),var(1,:)',var(2,:)',var(3,:),0.02,'nointerp');
% light;
var = velocity{ind};  
mag = sqrt(sum(var.^2,1));  cutoff = quantile(mag, 1);
var(:,mag > cutoff) = NaN;
qv = quiver3(Mesh(ind).verts(:,1),Mesh(ind).verts(:,2),Mesh(ind).verts(:,3),var(1,:)',var(2,:)',var(3,:)',0.75,'k','linewidth',2);
colorbar;  %caxis([0.3 0.6]); 
ylim([-7.5 -2.75]);
zlim([12 18]);
xlim([30 35]);

% export_fig(gcf,'C:\Users\rudi\Desktop\RD\pape main results\figures\power calc.png','-r400','-transparent')
        %%
        [~, A_force] = compute_force_integral(Mesh, matrix_props, input); % could make a very similar function that doesn't copy the values 3 times but who cares....
        A_scalar = squeeze(A_force(1,1:3:end,:));  % A_force has 3 copies of each entry for integrating fx, fy, fz components of force all the same way but here we are integrating a scalar
        if size(A_force,3) == 1
            A_scalar = A_scalar';
        end
        
        %         then use integrate unity similarly to how integrate force works - instead of A * f, it will be A * (f dot v)
        
        % basically,
        clear power
        for i_mesh = 1:size(A_scalar,2)
            power(i_mesh) = [A_scalar(:,i_mesh)]'  *  f_dot_v;  %cleverly, each page of A_scalar (going with each submesh) has zeros for all the verts not part of the current submesh, so they can each be multiplied by full f
        end
        %         eff = 6*pi*(  (Mesh(1).Volume * 3/4 / pi)^(1/3) ) *input.constants.mu * fits.converged.speed^2 / sum(power)
        powers(phase_ind) = sum(power);
        
      
    end
    
    switch input.kinematics_interp_method
        case 'spline'
            interpolant_power = spline([Solutions.phase(filter) ]',powers(filter)');
            power_mean = 1/diff(phase_bounds) * integral(@(x) fnval(interpolant_power, x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12);
        case 'trig'
            interpolant_power = trig_interp_fit(Solutions.phase,powers);
            power_mean = 1/diff(phase_bounds) * integral(@(x) trig_interp_eval(interpolant_power, x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12);
            
    end
    
    
    switch input.kinematics_interp_method
        case 'trig'
            mean_u = 1/diff(phase_bounds) * integral(@(x) trig_interp_eval(interpolant(end).U(1), x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12);
        case 'spline'
            mean_u = 1/diff(phase_bounds) * integral(@(x) fnval(interpolant(end).U(1), x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12);
    end
    
    eff.long_range = 6*pi*(  (Mesh(1).Volume * 3/4 / pi)^(1/3) ) *input.constants.mu * fits.converged.speed^2 / power_mean;
    eff.short_range = 6*pi*(  (Mesh(1).Volume * 3/4 / pi)^(1/3) ) *input.constants.mu * mean_u^2 / power_mean;
    
    
end

eff
%