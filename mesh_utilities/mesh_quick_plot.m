

bugtype = 'bacteria';
bugtype = 'dino';
AR1 = 9;  AR2 = 0.1;


folder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_3_modded\';
% folder = 'C:\Users\rudi\Desktop\RD\fucked\';
folder = 'E:\Hull\dinoflagellate\meshes_parallel2_3_modded\';

folder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_all\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_interp_4_initial\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_3\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_modded\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_interpBCs\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_static_tail_interped\';
% folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_orig\';
% folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_flipped_tail\';
% folder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_flipped2_interped\';
% folder = 'C:\Users\rudi\Desktop\RD\meshes_fat_long_sheet_tail\';
% folder = 'C:\Users\rudi\Desktop\RD\meshes_finer\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_finer_partial_turns2\';
% folder = 'C:\Users\rudi\Desktop\RD\centered tail\meshes_thin_tail - interp\';
folder = 'C:\Users\rudi\Desktop\RD\meshes_hairs_2_1.5_phases/';
folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_4_1\processed\';
folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_4_1\original\';
folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\45 deg flipped tail\processed\';
folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\regular\final\';
folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\unflipping tail\interpolated\';
% folder = 'C:\Users\rudi\Desktop\RD\opt_meshes\';
% folder = 'C:\Users\rudi\Desktop\RD\meshes_1x_thin_tail_length\';
% folder = 'C:\Users\rudi\Desktop\RD\phases_hairs_3.5_1.5\orig\';
% folder = 'C:\Users\rudi\Desktop\RD\phases_hairs_3.5_1.5\computed BCs\';
folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\stopped transverse 45 deg flipped tail\';
folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\regular\cleaned\';
folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\centered_tail\interpolated\';
folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\final\';

folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\6x tail length\original\';
% folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\stopped transverse\';

actively_flipping_tail = false;
n_phase_pts = 128;


switch bugtype
    case 'dino'
        
        if actively_flipping_tail
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

        
        step = 2;  % 109 110 112 113 115
        % 30 fucked
        
%         submeshes = [1:6];
%           submeshes = [ 1 3];
        %    submeshes =[ 1 2 4:6];
        %    submeshes = [2 4 5 6];
%           submeshes = [2 4:6];
         submeshes = [1 3];
        
        
        time = time_list(step_list == step);
        phase = phase_list(step_list == step);
        basename = ['_step_',num2str(step),'_time_', sprintf('%0.15f',time),'_phase_',sprintf('%0.15f',phase)];
        
        
        % names = {'Transverse','Wingtip'};
        names = {'Body', 'Transverse','Tail', 'Coplanar_Hairs' , 'Normal_Top_Hairs'  , 'Normal_Bottom_Hairs'};
        
    case 'bacteria'
        names = {''};
        basename = ['curved_rod_' ,'AR1_',num2str(AR1),'_AR2_',num2str(AR2)] ;
        submeshes = 1;
end

clear shat Mesh Metadata
for i = 1:length(names)
    
    
    
    %     [Mesh(i), Metadata] = load_mesh(['E:\Hull\dinoflagellate\meshes_fat_tail\',names{i},'_step_7_time_0.001188858695652_phase_0.343611696486384.dat'],[],[],'mesh');
    %    [Mesh(i), Metadata] = load_mesh(['E:\Hull\dinoflagellate\meshes_transverse_hairs_wider_parallel\',names{i},'_step_0_time_0.000000000000000_phase_0.000000000000000.dat'],[],[],'mesh');
    
    %     [Mesh(i), Metadata] = load_mesh(['E:\Hull\dinoflagellate\meshes_transverse_hairs_wide\',names{i},'_step_0_time_0.000000000000000_phase_0.000000000000000.dat'],[],[],'mesh');
    % [Mesh(i), Metadata] = load_mesh(['E:\Hull\dinoflagellate\meshes_transverse_hairs_parallel2_3_all\',names{i},'_step_0_time_0.000000000000000_phase_0.000000000000000.dat'],[],[],'mesh');
    
    %[Mesh(i), Metadata] = load_mesh(['E:\Hull\dinoflagellate\meshes_parallel2_3_all\',names{i},'_step_4_time_0.000679347826087_phase_0.196349540849362.dat'],[],[],'mesh');
    [Mesh(i), Metadata(i)] = load_mesh([folder,names{i},basename,'.dat'],[],[],'mesh');
    
    % [Mesh(i), Metadata] = load_mesh(['E:\Hull\dinoflagellate\meshes_parallel2_3\',names{i},'_step_70_time_0.011888586956522_phase_3.436116964863836.dat'],[],[],'mesh');
    
    
    %    [Mesh(i), Metadata] = load_mesh([mesh_file],[],[],'mesh');
 
end

for i = 1:length(names)
Mesh(i).name = names{i};
end
%%

%   mesh_file = 'E:\Hull\meshes\ellipsoid_AR1_0.6_AR2_7.dat';
%
%   mesh_file = 'E:\Hull\all_meshes\curved_rod_AR1_5_AR2_0.65.dat';
%
% %   mesh_file = 'C:\Users\rudi\Desktop\RD\swept_meshes\curved_rod_AR1_7.5_AR2_0.75.dat';
% %
% %   mesh_file = 'E:\Hull\dinoflagellate\meshes_transverse_hairs_wide\Transverse_step_0_time_0.000000000000000_phase_0.000000000000000.dat';
% %
%     mesh_file = 'E:\Hull\Oscar_meshes\tail_radius_0.012_amp_0.2_lambda_2.34_nlambda_2.564102564102564.dat';
%
%     mesh_file = 'C:\Users\rudi\Desktop\RD\meshes_parallel2_uberinterp\
%
%  i = 1;  clear Mesh Metadata
% % %  [Mesh(i), Metadata] = load_mesh([mesh_file],[],[],'mesh');
%    [Mesh(i), Metadata] = load_mesh([mesh_file],[],[]);
%764  3424  1738

[Mesh] = global_inds(Mesh);
[Mesh] = renumber_Mesh(Mesh);

clear temp_input
temp_input.performance.nthreads = 8;
temp_input.accuracy.integration_tol.area.abstol = 1E-6;
temp_input.accuracy.integration_tol.area.reltol = 1000;
temp_input.accuracy.integration_tol.area.maxevals = Inf;

temp_input.accuracy.integration_tol.centroid.abstol = 1E-6;
temp_input.accuracy.integration_tol.centroid.reltol = 1000;
temp_input.accuracy.integration_tol.centroid.maxevals = Inf;

temp_input.accuracy.integration_tol.volume.abstol = 1E-6;
temp_input.accuracy.integration_tol.volume.reltol = 1000;
temp_input.accuracy.integration_tol.volume.maxevals = Inf;

[Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);

Mesh(3)  =  shiftMesh(Mesh(3),geom.tail.shift');  

%%
% figure;  plot_mesh(Mesh(submeshes),[3 3 3 3 3 3]);   axis tight; % light;
figure(725); clf;  [s,e] = plot_mesh(Mesh(submeshes),[2 2 2 4 2 2]);   axis tight; % light;
% set(gca,'view',[0 0]); % side
set(gca,'view',[90 0]); % front
title(folder,'interpreter','none');
set(e,'edgealpha',0.1);
   
    light
% figure(586); clf;  plot_mesh(Mesh(4),[3]);   axis tight;  light;
% fucked = 11 17 34 35 37 38 39 40 54 55 60 61 82 83 98 99 100
%%
  export_fig(gcf,['C:\Users\rudi\Desktop\RD\pape main results\figures\','Body 6x Tail','.png'],'-r400','-transparent');

% figure;  plot_mesh(Mesh(2));   light