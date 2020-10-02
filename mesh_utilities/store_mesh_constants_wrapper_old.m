function Mesh = store_mesh_constants_wrapper(Mesh, input)


%need to deal with incredibly annoying mex rules by removing possible new
%fields and adding back later, in case store_mesh_constants is ever done after
%submeshes have been combined, etc
expected = {'n_vert','n_elem','verts','elems','indices'};
fields = fieldnames(Mesh);
extra = setdiff(fields,expected); %store_mesh_constants_mex isn't going to expect these fields and will fail if it sees them
for f = 1:length(extra)
    for j = 1:length(Mesh)
        temp(j).(extra{f}) = Mesh(j).(extra{f});
    end
end
    
Mesh = rmfield(Mesh,extra);
    
[Mesh] = store_mesh_constants_init(Mesh);

% [Mesh] = store_mesh_constants_mexed(Mesh,input);
[Mesh] = store_mesh_constants_mexed(Mesh,input);

%add extra fields back before returning
for f = 1:length(extra)
    for j = 1:length(Mesh)
        Mesh(j).(extra{f}) = temp(j).(extra{f});
    end
end

