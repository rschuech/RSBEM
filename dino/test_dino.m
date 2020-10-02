clear Mesh
cd C:\Users\rudi\Desktop\RD\dino_code\
addpath(genpath('./')); %add all subfolders to path so that subfunctions will be found

input.paths.datfolder = 'C:\Users\rudi\Desktop\RD\meshes_grouped_fast\';
input.paths.datfolder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail_nosphere2\';
%input.paths.datfolder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail_nosphere_big_groove\';
input.paths.datfolder = 'C:\Users\rudi\Desktop\RD\meshes_thin_tail_nosphere_tail_angle_25\';
input.paths.datfolder = 'C:\Users\rudi\Desktop\RD\meshes_thick_tail_nosphere_tail_angle_0\';
input.paths.datfolder = 'E:\Hull\dinoflagellate\meshes_mark5\';

input.paths.namebase.body = 'Body_step_0_time_0.0_phase_0.0';
input.paths.namebase.transverse = 'Transverse_step_0_time_0.0_phase_0.0';
input.paths.namebase.tail = 'Tail_step_0_time_0.0_phase_0.0';
input.bugtype = 'dino';
input.problemtype = 'freeswim';
input.performance.nthreads = 8;
input.tail.motorBC = char(zeros(1,0));
input.ignore_interaction = false;
input.accuracy.eps2 = (4E-5  ) ^2;   %tail is 5 times bigger than orig, although body is 15 times bigger now, need to keep tail integrals accurate
input.constants.mu = 1E-3 / 1E6; %viscosity in kg / (s micron)
%gets multiplied by boundary integral
input.constants.multfactor = 1/8/pi;
input.performance.debug_mode = false;
input.output.interpolation.fignum = 123;
input.accuracy.interpolation.max_iter = 6;
input.output.interpolation.doplot = true;
input.output.interpolation.doprint = true;
input.paths.dumpfolder = 'C:\Users\rudi\Desktop\RD\dino_dumps\';
input.paths.fullnamebase = 'correct_scale2';
input.accuracy.interpolation.max_rel_diff = 0.1;
%%

%time = 6.28318530718;
time = 0;

load([input.paths.datfolder,'Metadata_step_0_time_0.0_phase_0.0.mat']);

for bi = 1
    clear Mesh
    switch bi
        
        case 1 % body + tail + transverse
            str = 'body_tail_transverse';
            Mesh(1) = load_mesh([input.paths.datfolder,input.paths.namebase.body,'.dat'],[],[Metadata.body.rand_inds],'mesh');
            Mesh(2) = load_mesh([input.paths.datfolder,input.paths.namebase.transverse,'.dat'],[],[Metadata.transverse.rand_inds],'mesh');
            Mesh(3) = load_mesh([input.paths.datfolder,input.paths.namebase.tail,'.dat'],[],[Metadata.tail.rand_inds],'mesh');
            
            Mesh(1).name = 'body';
            Mesh(2).name = 'transverse';
            Mesh(3).name = 'tail';
            
        case 2  %body + tail
            str = 'body_tail';
            
            Mesh(1) = load_mesh([input.paths.datfolder,input.paths.namebase.body,'.dat'],[],[Metadata.body.rand_inds],'mesh');
            Mesh(2) = load_mesh([input.paths.datfolder,input.paths.namebase.tail,'.dat'],[],[Metadata.tail.rand_inds],'mesh');
            
            Mesh(1).name = 'body';
            Mesh(2).name = 'tail';
        case 3 %body + transverse
            str = 'body_transverse';
            Mesh(1) = load_mesh([input.paths.datfolder,input.paths.namebase.body,'.dat'],[],[Metadata.body.rand_inds],'mesh');
            Mesh(2) = load_mesh([input.paths.datfolder,input.paths.namebase.transverse,'.dat'],[],[Metadata.transverse.rand_inds],'mesh');
            
            Mesh(1).name = 'body';
            Mesh(2).name = 'transverse';
        case 4 %tail + transverse
            str = 'tail_transverse';
            Mesh(1) = load_mesh([input.paths.datfolder,input.paths.namebase.transverse,'.dat'],[],[Metadata.transverse.rand_inds],'mesh');
            Mesh(2) = load_mesh([input.paths.datfolder,input.paths.namebase.tail,'.dat'],[],[Metadata.tail.rand_inds],'mesh');
            
            Mesh(1).name = 'transverse';
            Mesh(2).name = 'tail';
        case 5 %tail
            str = 'tail';
            Mesh(1) = load_mesh([input.paths.datfolder,input.paths.namebase.tail,'.dat'],[],[Metadata.tail.rand_inds],'mesh');
            Mesh(1).name = 'tail';
        case 6 %transverse
            str = 'transverse';
            Mesh(1) = load_mesh([input.paths.datfolder,input.paths.namebase.transverse,'.dat'],[],[Metadata.transverse.rand_inds],'mesh');
            Mesh(1).name = 'transverse';
    end
    
    
    str
    
    for si = 1:length(Mesh)
        Mesh(si).orientation = [1 0 0; 0 1 0; 0 0 1;]';
        Mesh(si).refpoints = [0 0 0]';
    end
    
    input.accuracy.integration_tol.traction.abstol =   16E-3  *4  ; %scale by min (worst case) element size compared to convergence test case
    input.accuracy.integration_tol.traction.reltol = 0;  %control error via abstol
    input.accuracy.integration_tol.traction.maxevals = Inf; %don't limit maxevals
    input.accuracy.integration_tol.force.abstol = 1E-4  *4 ; %1E-7
    input.accuracy.integration_tol.force.reltol = 0;
    input.accuracy.integration_tol.force.maxevals = Inf;
    input.accuracy.integration_tol.torque.abstol = 1E-4  *4; %1E-5
    input.accuracy.integration_tol.torque.reltol = 0;  %2.5E-5
    input.accuracy.integration_tol.torque.maxevals = Inf;
    input.accuracy.integration_tol.area.abstol = 0 ;
    input.accuracy.integration_tol.area.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
    input.accuracy.integration_tol.area.maxevals = Inf;
    input.accuracy.integration_tol.volume.abstol = 0 ;
    input.accuracy.integration_tol.volume.reltol = 1E-1  ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
    input.accuracy.integration_tol.volume.maxevals = Inf;
    input.accuracy.integration_tol.centroid.abstol = 0 ;
    input.accuracy.integration_tol.centroid.reltol = 10 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
    input.accuracy.integration_tol.centroid.maxevals = Inf;
    
    
    
    [Mesh] = global_inds(Mesh);
    [Mesh] = renumber_Mesh(Mesh);
    % create smaller input struct with just necessary parameters to reduce mex
    % pain
    temp_input.performance.nthreads = 20;
    
    
    temp_input.accuracy.integration_tol.area.abstol = input.accuracy.integration_tol.area.abstol ;
    temp_input.accuracy.integration_tol.area.reltol = input.accuracy.integration_tol.area.reltol ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
    temp_input.accuracy.integration_tol.area.maxevals = input.accuracy.integration_tol.area.maxevals;
    temp_input.accuracy.integration_tol.centroid.abstol = input.accuracy.integration_tol.centroid.abstol ;
    temp_input.accuracy.integration_tol.centroid.reltol = input.accuracy.integration_tol.centroid.reltol; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
    temp_input.accuracy.integration_tol.centroid.maxevals = input.accuracy.integration_tol.centroid.maxevals;
    temp_input.accuracy.integration_tol.volume.abstol = input.accuracy.integration_tol.volume.abstol ;
    temp_input.accuracy.integration_tol.volume.reltol = input.accuracy.integration_tol.volume.reltol ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
    temp_input.accuracy.integration_tol.volume.maxevals = input.accuracy.integration_tol.volume.maxevals;
    
    [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
    
    
    
    
    clear matrix_props
    % matrix_props.n_col = sum([Mesh.n_vert]);  %number collocation points in A matrix, could be greater than number of points at which to determine traction, in theory
    %    matrix_props.n_col = Mesh(end).indices.glob.vert( find( ~isnan(Mesh(end).indices.glob.vert), 1, 'last' ) );  %highest global index must be in last submesh, last nonNaN entry since it must monotonically increase
    matrix_props.n_col = Mesh(end).indices.glob.unq_bounds.vert(2);
    
    if strcmp(input.problemtype,'freeswim')
        switch input.bugtype
            case 'bacteria'
                switch input.tail.motorBC
                    case 'freq'
                        matrix_props.n_rows = matrix_props.n_col * 3 + 6;  %usual unknowns plus 3 translation components and 3 rotation components of body
                    case 'torque'
                        matrix_props.n_rows = matrix_props.n_col * 3 + 7;  %usual unknowns plus 3 translation components and 3 rotation components of body + rotation rate of tail
                end
            case 'dino'
                matrix_props.n_rows = matrix_props.n_col * 3 + 6;  %usual unknowns plus 3 translation components and 3 rotation components of body
        end
    else
        matrix_props.n_rows = matrix_props.n_col * 3;
    end
    matrix_props.n_cols = matrix_props.n_rows; %number of columns in A matrix
    
    matrix_props.Col_inds = save_Col_inds_mexed(Mesh,input.performance.nthreads);
    
    %%
    
    clear assembly_input
    assembly_input.performance.nthreads = input.performance.nthreads;
    assembly_input.problemtype = input.problemtype;
    assembly_input.tail.motorBC = input.tail.motorBC;
    assembly_input.ignore_interaction = input.ignore_interaction;
    assembly_input.accuracy.integration_tol.traction = input.accuracy.integration_tol.traction;
    assembly_input.accuracy.integration_tol.force = input.accuracy.integration_tol.force;
    assembly_input.accuracy.integration_tol.torque = input.accuracy.integration_tol.torque;
    assembly_input.accuracy.eps2 = input.accuracy.eps2;
    assembly_input.constants.multfactor = input.constants.multfactor;
    assembly_input.constants.mu = input.constants.mu;  %not actually needed for matrix assembly, but needed for yprime() during timestepping interpolation
    assembly_input.performance.debug_mode = input.performance.debug_mode;
    assembly_input.skip_rigid_integrals = false;
    assembly_input.skip_traction_integrals = false;
    assembly_input.bugtype = 'dino';
    
    ass_tic = tic;
    [ A, ~,A_force,A_torque,~] = matrix_assembly_mex_wrapper( Mesh,matrix_props,assembly_input );
    total_assembly = toc(ass_tic)
    
        Answers_thick_tail_angle_0.(str).A_force = A_force;
    Answers_thick_tail_angle_0.(str).A_torque = A_torque;
  
    
    %%
    % flowcase = 'rx';
    % U0 = 50;  %microns / s    arbitrary, for calculation of friction coeffs
    % Omega0 = 1;  %rad / s     arbitrary, for calculation of rotational friction coeffs
    % BCs = set_forced_BCs(U0, Omega0,flowcase); %defines BCs struct with U and Omega vectors
    % input.problemtype = 'forced';
    % [f] = matrix_solve(Mesh, input, matrix_props, A(1:end-6,1:end-6), BCs); %solves matrix equation and returns solution vector, in this case just traction
    % [forces] = integrate_traction(A_force,A_torque,f);  %integrates traction and returns forces and torques on entire object
    % [fcoeffs_temp] = friction_coeffs(forces,BCs);  %computes translational and rotational friction coeffs
    
    for n = 1:length(Mesh)
        name = Mesh(n).name;
        [~,inds] = ismember(Mesh(n).indices.orig.vert,Metadata.(name).indices.orig.vert);  %both current Mesh and Mesh used in Metadata were randomized so need to reorder BCs to match current order
        BCs.freeswim.(name) = Metadata.(name).BCs(inds,:);
    end
    
    solve_tic = tic;
    [RHS] = assemble_RHS(input,Mesh, matrix_props, BCs);
    %[f] = matrix_solve(input, matrix_props, A(1:end-6,1:end-6), RHS); %solves matrix equation and returns solution vector, in this case just traction
    [f, kinematics] = matrix_solve(input, matrix_props, A, RHS); %solves matrix equation and returns solution vector, in this case just traction
    solve_time = toc(solve_tic)
    
     Answers_thick_tail_angle_0.(str).f = f;
    Answers_thick_tail_angle_0.(str).kinematics = kinematics;

    
    
    
end
stopooo


[forces] = integrate_traction(A_force,A_torque,f);  %integrates traction and returns forces and torques on entire object
[fcoeffs_temp] = friction_coeffs(forces,BCs);  %computes translational and rotational friction coeffs

%% copied from matrix_solves
vertlist = [];
for i = 1:length(Mesh)
    inds = find( Mesh(i).indices.glob.vert >= Mesh(i).indices.glob.unq_bounds.vert(1) &   Mesh(i).indices.glob.vert <= Mesh(i).indices.glob.unq_bounds.vert(2)); %local inds of "global" verts belonging to this mesh
    vertlist = [vertlist; Mesh(i).verts(inds,:)];
end
[u] = compute_BCs_mexed(BCs,vertlist,input.performance.nthreads);
%%
u =  RHS(1:matrix_props.n_col*3);
%%
f2 = f;  %only take part of f for free swimming since last few are not forces
f2m = reshape(f2,3,[]);
f2mag = sqrt(sum(f2m.^2));
um = reshape(u,3,[]);
umag = sqrt(sum(um.^2));

delta = 1;

var = 'f2';  %either u or f2, speed or force
eval(['var1 = ',var,'m',';']);
eval(['var2 = ',var,'mag',';']);

n_refines = 2;
Surface_vis = refine_vis_surface(Mesh,n_refines);
clear var_vis* varM
for i = 1:length(Mesh)
    Surface_vis(i).phi = NaN(Surface_vis(i).n_vert,6);
    for vi = 1:Surface_vis(i).n_vert
        [~, ~, Surface_vis(i).phi(vi,:), ~] = T6interp(NaN(6,3),  Surface_vis(i).xi_eta(vi,1),  Surface_vis(i).xi_eta(vi,2),  [Mesh(i).elem_params(1,Surface_vis(i).xi_eta_elem(vi)),  Mesh(i).elem_params(2,Surface_vis(i).xi_eta_elem(vi)),  Mesh(i).elem_params(3,Surface_vis(i).xi_eta_elem(vi))]   );
    end
    
    varM{i} = var1(:,Mesh(i).indices.glob.vert);  %f2m needs to be copied into local verts for each submesh since they share some verts and thus f values
    
    var_vis{i} = NaN(3,Surface_vis(i).n_vert);
    for vi = 1:Surface_vis(i).n_vert
        orig_vals = varM{i}(:, Mesh(i).elems(Surface_vis(i).xi_eta_elem(vi),:));  %f at verts of orig element
        var_vis{i}(:,vi)   =  orig_vals * Surface_vis(i).phi(vi,:)';
    end
    var_vis_mag{i} = sqrt(sum(var_vis{i}.^2))';
    
end

%p = patch('faces',Mesh2.elems,'vertices',Mesh2.verts,'facevertexcdata',[colors1(colorinds1,:); colors2(colorinds2,:)],'facecolor','interp','edgecolor','none','facealpha',facealpha);

facealpha = 1;

sidepts = 7;  %number of points along triangle edges was 7

Edge_vis = refine_vis_edges(Mesh,sidepts);
%%
clf
clear colors colorinds
mincolor = min(vertcat(var_vis_mag{:})) * 1;
maxcolor = max(vertcat(var_vis_mag{:})) * 1;

%    [colors{1},colorinds{1}] = colordata(500,'jet',[mincolor maxcolor],var_vis_mag{1});
%         [colors{2},colorinds{2}] = colordata(500,'jet',[mincolor maxcolor],var_vis_mag{2});
%         [colors{3},colorinds{3}] = colordata(500,'jet',[mincolor maxcolor],var_vis_mag{3});

[colors{1},colorinds{1}] = colordata(500,'jet',[min(var_vis_mag{1}) max(var_vis_mag{1})*0.1],var_vis_mag{1});
[colors{2},colorinds{2}] = colordata(500,'jet',[min(var_vis_mag{2}) max(var_vis_mag{2})*0.5],var_vis_mag{2});
[colors{3},colorinds{3}] = colordata(500,'jet',[min(var_vis_mag{3}) max(var_vis_mag{3})*0.8],var_vis_mag{3});
edgealpha = 0.15;
clear s
for i = 1:length(Surface_vis)
    s(i) = patch('faces',Surface_vis(i).elems,'vertices',Surface_vis(i).verts,'facevertexcdata',colors{i}(colorinds{i},:),'edgecolor','none','facecolor','interp','facealpha',facealpha);
    hold on
    % q(i) = quiver3(Surface_vis(i).verts(:,1),Surface_vis(i).verts(:,2),Surface_vis(i).verts(:,3),var_vis{i}(1,:)',var_vis{i}(2,:)',var_vis{i}(3,:)',1,'k');
    e(i) = patch('faces',Edge_vis(i).elems,'vertices',Edge_vis(i).verts,'facecolor','none','edgecolor','k','edgealpha',edgealpha);
end
hold off
axis equal; grid on;  box on;
xlabel('x');  ylabel('y');  zlabel('z');
set(s,'facealpha',1)
light








%% interp stepping
folder = 'C:\Users\rudi\Desktop\RD\meshes_grouped_fast\';

clear Inds Files Times Phases
names = {'Body','Transverse','Tail','Metadata'};
for n = 1:length(names)
    name = names{n};
    
    files = dir([folder,name,'*']);
    files = {files.name};
    %files = files(3:end);
    
    clear times phases
    for i = 1:length(files)
        ind1 = strfind(files{i},'time_') + 5;
        ind2 = strfind(files{i},'_phase') - 1;
        ind3 = strfind(files{i},'.dat') - 1;
        times(i) = str2double(files{i}(ind1:ind2));
        phases(i) = str2double(files{i}(ind2+8:ind3));
    end
    
    [~,inds] = sort(times);
    Inds{n} = inds;
    Files{n} = files(inds);
    Times{n} = times(inds);
    Phases{n} = phases(inds);
end

if ~isequal(Times{1},Times{2},Times{3})
    disp('problemo')
    pause
else
    Time = Times{1};
    Phase = Phases{1};
end

Mesh_files.time = Time;
Mesh_files.phase = Phase;
Mesh_files.body = Files{1};
Mesh_files.transverse = Files{2};
Mesh_files.tail = Files{3};
Mesh_files.metadata = Files{4};
%%


BCs.freeswim.phase_speed = 1.5;  % lowest common period T = 4*pi s, covered in 2*pi rad, is 1/2 rad/sec for the entire beat cycle
y0 = [ [NaN NaN NaN]' ; [NaN NaN NaN ]' ]; %actually isn't used for dino
y0 = [Mesh(1).refpoints(:,1)' 0 0 0 ];
recycle_A = false;
interp_y_tic = tic;
interp_y;
timings.interpolation = toc(interp_y_tic);


%% unroll all fields and vector components of last and most accurate interpolant
cc = 0;  clear best_interpolant
fields = fieldnames(interpolant(end));
for f = 1:length(fields)  %each vector variable, e.g. U, Omega, omega
    for i = 1:length( interpolant(end).(fields{f}) )  %components of each variable, e.g U(1:3), Omega(1:3)
        cc = cc+1;
        best_interpolant(cc) = interpolant(end).(fields{f})(i);
    end
end
