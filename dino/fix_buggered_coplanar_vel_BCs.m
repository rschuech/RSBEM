


folder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_3_1.5\final\';

files = dir([folder,'Metadata*.mat']);  files = {files.name};

parfor f = 1:length(files)
%     temp = load([folder,files{f}]);
%     Metadata = temp.Metadata;
    
    temp = matfile([folder,files{f}],'writable',true);
    Metadata = temp.Metadata;
    [~, vel ,~,~] = transverse_hairs_parameterized(Metadata.Coplanar_Hairs.parameters(:,1),Metadata.Coplanar_Hairs.parameters(:,2), Metadata.time, Metadata.geom.transverse , 'Coplanar_Hairs' );
    
    
    Metadata.Coplanar_Hairs.BCs = vel';
    temp.Metadata = Metadata;
    
%     save([folder,files{f}],'Metadata');
    
end