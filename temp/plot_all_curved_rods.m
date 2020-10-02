
folder = 'C:\Users\rudi\Desktop\RD\opt_meshes\';

for i = 1:length(Results)
mesh_file = ['C:\Users\rudi\Desktop\RD\swept_meshes\curved_rod_AR1_',num2str(Results(i).AR1),'_AR2_',num2str(Results(i).AR2),'.dat'] 
submeshes = 1;

clear shat Mesh Metadata

%%

% % %  [Mesh(i), Metadata] = load_mesh([mesh_file],[],[],'mesh');
   [Mesh, Metadata] = load_mesh([mesh_file],[],[]);

 
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

Volumes(i) = Mesh.Volume;
AR1_AR2_Results(i,:) = [Results(i).AR1 Results(i).AR2];
AR1_AR2_Metadata(i,:) = [Metadata.geom.AR1 Metadata.geom.AR2];

continue

figure(46); 
clf
plot_mesh(Mesh(submeshes),[3 3 3 3 3 3]);   axis tight;  light;
% set(gca,'view',[0 0]);

pause

end
% figure(586); clf;  plot_mesh(Mesh(4),[3]);   axis tight;  light;
% fucked = 11 17 34 35 37 38 39 40 54 55 60 61 82 83 98 99 100
%%

% figure;  plot_mesh(Mesh(2));   light