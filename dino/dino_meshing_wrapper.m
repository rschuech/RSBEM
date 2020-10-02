runtype = 'meshes';  % 'meshes' for original Salome mesh generation, 'BCs' for compute_BCs script that generates velocity BCs and adjusts mesh verts

tail_mesh_only = true; % skip hair sheets geom and meshes and slow Partition operation in Salome, only output tail mesh
static_transverse = false;

max_time = 600;  % how long to allow Salome to run for (sec); 600 seems OK
% max_time = 20;
check_freq = 10; % how frequently to check whether mesh files have been created (sec)

num_parallel_sessions = 20; % how many Salomes to run in parallel at any moment

coplanar_length = 2;  parallel_hair_length_input = coplanar_length;  % 3.5
normal_length = 1.5;  normal_hair_length_top_input = normal_length;  normal_hair_length_bot_input = normal_length;  % 1.5

parallel_hairs_nsegs_radial_input = 11;  %15
normal_hairs_nsegs_radial_input = 7;  % 4

% coplanar_length normal_length coplanar_nsegs normal_nsegs
% 4               1               25             4
% 2             1.5               11             7
% 2.5           1.5                9             7
% 3             1.5               11             7

n_phase_pts = 32;

% phase_indices = (0:127);
phase_indices = setdiff(0:31, []);
% phase_indices = 1;
%  phase_indices = [0 1 2 3 63 ];
%  phase_indices = 0;
 
%  folder = ['phases_hairs_',num2str(coplanar_length),'_',num2str(normal_length),  '/','unflipping tail/'];
%  folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2.25_1.5\';
  folder = ['phases_hairs_',num2str(coplanar_length),'_',num2str(normal_length),  '/'];
  folder = ['6x tail length/'];
 
    digits = 16;
    

    paths = [];
    paths.template_file = 'dino_template.py';
 
     switch getenv('computername')
        case 'UBERTOP'
            run_bat_template = 'C:\Hull\dinoflagellate\mesh_automation\run_salome_template.bat';
            run_py_template = 'C:\Hull\dinoflagellate\mesh_automation\runSalome_template.py';
            kill_salome_carefully_template = 'C:\Hull\dinoflagellate\mesh_automation\kill_salome_carefully_template.py';
            paths.outfolder = 'C:/Hull/dinoflagellate/mesh_automation/';  %be sure to use / instead of \ here; regexprep won't work with \
            
            paths.infolder = 'C:/Hull/dinoflagellate/mesh_automation/';
        case {'CFD01','CFD02','CFD03','CFD04'}
            run_bat_template = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\run_salome_template.bat';
            run_py_template = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\runSalome_template.py';
%             kill_salome_carefully_template = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\kill_salome_carefully_template.py';
            
            paths.outfolder = 'C:/Users/rudi/Desktop/RD/dino_mesh_automation/';
            paths.infolder = 'C:/Users/rudi/Desktop/RD/dino_mesh_automation/';
            
%             paths.outfolder = 'C:\Users\rudi\Desktop\RD\recessed hairs\';
% paths.outfolder = folder;
%             paths.infolder = paths.outfolder;
%             folder_original = folder;
%             folder_processed = folder;
    end
    
    outputprefix = 'dino';  %put this at beginning of output file names
    suffix = '';  %append at end of file names, if used, start with _
    
    folder_input = [paths.outfolder];
    if ~exist(folder_input,'dir')
        mkdir(folder_input);
    end
        folder_temp = [paths.outfolder, folder,'temp/'];
    if ~exist(folder_temp,'dir')
        mkdir(folder_temp);
    end
    folder_original = [paths.outfolder,folder,'original/'];
    if ~exist(folder_original,'dir')
        mkdir(folder_original);
    end
    folder_processed = [paths.outfolder,folder,'processed/'];
    if ~exist(folder_processed,'dir')
        mkdir(folder_processed);
    end
     folder_interpolated = [paths.outfolder,folder,'interpolated/'];
    if ~exist(folder_interpolated,'dir')
        mkdir(folder_interpolated);
    end
    
    
    switch runtype
        case 'meshes'
            %%
            !taskkill /FI "IMAGENAME eq SALOME*" /T /F
            !taskkill /FI "IMAGENAME eq omniNames*" /T /F
            temp = dir('C:\Users\rudi\AppData\Local\Temp\**\.salome_PortManager*');
            for tp = 1:length(temp)
                delete([temp(tp).folder,'/',temp(tp).name]);
            end
             [gar,bage] = system('C:\Users\rudi\Desktop\RD\Salome\SALOME-8.3.0-WIN64\WORK\kill_salome.bat');
             %%
            gen_dino_meshes;
            
        case 'BCs'
            while true
                compute_dino_BCs;
                pause(10);
            end
    end
    
    