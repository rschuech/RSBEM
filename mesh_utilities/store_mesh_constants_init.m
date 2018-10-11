function [Mesh] = store_mesh_constants_init(Mesh)


for j = 1:length(Mesh)
%     for i = 1:N
%         Mesh2(j).(fields{i}) = Mesh(j).(fields{i});
%     end
    Mesh(j).area = NaN(size(Mesh(j).elems,1),1);
    Mesh(j).Centroid = NaN(3,1);
    Mesh(j).elem_params = NaN(3,size(Mesh(j).elems,1));
    Mesh(j).centroid = NaN(3,size(Mesh(j).elems,1));
    Mesh(j).Volume = NaN;
end