folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\stopped transverse\temp\';
step_0_folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\stopped transverse\temp\step 0\';
outfolder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\stopped transverse\';

template = '_step_0_time_0.000000000000000_phase_0.000000000000000';

names = {'Body','Coplanar_Hairs','Normal_Bottom_Hairs','Normal_Top_Hairs','Tail','Transverse','Metadata'};


for n = 1:length(names)
    name = names{n};
    
    files = dir([folder,name,'*']);  files = {files.name};
    
    
    for i = 1:length(files)
        
        switch name
            case 'Metadata'
                suffix = '.mat';
            otherwise
                suffix = '.dat';
        end
        copyfile([step_0_folder,name,template,suffix] , [outfolder,files{i}]);
        
        
    end
    
    
end

%%
n_phase_pts = 128;  temp = 1;
step_list = 0:(n_phase_pts - temp);

total_period = 1 / 46.0 ;
time_list = linspace(0, total_period, n_phase_pts + 1)  ;
time_list = time_list(1:end - temp);  % leave off last value which is same as first
phase_list = linspace(0, 2*pi, n_phase_pts + 1);   % *always* 2 pi rad in a cycle
phase_list = phase_list(1:end - temp);  % leave off last value which is same as first

files = dir([outfolder,'Metadata','*']);  files = {files.name};

names = {'Body','Coplanar_Hairs','Normal_Bottom_Hairs','Normal_Top_Hairs','Transverse'};

for i = 1:length(files)
    file = files{i};
    ind = strfind(file,'step_');
    temp = sscanf(files{i}(ind+5:end),'%f_time_%f_phase_%f.dat'); %[step, time, phase]
    time = temp(2);
    index = temp(1);
    phase = temp(3);
    
    load([outfolder,file]);
    
    Metadata.time = time_list(step_list == index);
    
    for b = 1:length(names)
        Metadata.(names{b}).BCs(:) = 0;
    end
    
     save([outfolder,file],'Metadata');
    
end

%%


folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\stopped transverse 45 deg flipped tail\';
parameters_file = [folder,'Parameters_step_1_time_0.000169836956522_phase_0.049087385212341.txt'];
fat_tail = false;  coplanar_length = 2;  normal_length = 1.5;
dino_geom_parameters_from_file;

files = dir([folder,'Tail','*']);  files = {files.name};

%    tail_vel(vert_i,:) = rotate_arbitrary_axis(tail_vel(vert_i,:),[0 0 0]', [0 1 0]', pi/4);  %rotate around origin since this is a velocity vector
        
for i = 1:length(files)
    file = files{i};

    temp = dlmread([folder,file]);
clear Mesh
Mesh.n_vert = temp(1,1);  %# vertices, may or may not be the number that's actually used - see below
Mesh.n_elem = temp(1,2); %# elements

Mesh.indices.orig.vert = temp(2:Mesh.n_vert+1,1);
Mesh.indices.orig.elem = temp(Mesh.n_vert+2:end,1);

Mesh.verts = temp(2:Mesh.n_vert+1,2:4); %vertex coords

Mesh.elems = temp(Mesh.n_vert+2:end,3:end);  %vertex indices that form each element

 Mesh.verts = rotate_arbitrary_axis(Mesh.verts, geom.tail.translation, [0 1 0]', pi/4);
       
  delete([folder,file]);
 dlmwrite([folder,file],[Mesh.n_vert Mesh.n_elem],'-append','delimiter',' ','precision',16)
        dlmwrite([folder,file],[Mesh.indices.orig.vert    Mesh.verts],'-append','delimiter',' ','precision',16)
        dlmwrite([folder,file],[Mesh.indices.orig.elem      repmat(206,Mesh.n_elem,1)     Mesh.elems],'-append','delimiter',' ','precision',16)
       
        
end


%%
files = dir([folder,'Metadata','*']);  files = {files.name};
for i = 1:length(files)
    file = files{i};
clear Metadata
load([folder,file]);

   Metadata.Tail.BCs = rotate_arbitrary_axis(Metadata.Tail.BCs, [0 0 0]', [0 1 0]', pi/4);  %rotate around origin since this is a velocity vector
     save([folder,file],'Metadata');
end