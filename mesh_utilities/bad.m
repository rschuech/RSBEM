metaname = [meshname(1:end-4),'_metadata.mat']; %metadata filename should always correspond to mesh filename

if exist(metaname,'file')
    load(metaname);  %contains one variable, Metadata
    if isfield(Metadata.mesh,'topology') && ~strcmp(Metadata.mesh.topology, 'sphere')
        disp(['Mesh topologically equivalent to ',Metadata.mesh.topology]);
    end
else
    disp('Metafile doesn''t exist; aborting');
    Mesh = [];
    Metadata = [];
    return
end


   Metadata.mesh.meshing_succeeded = false;
            save(paths(c).metafile,'Metadata');
            succeeded(c) = false;
            reason(c) = -1;
            save(paths(path_ind).metafile,'Metadata');
            
            
            
            c