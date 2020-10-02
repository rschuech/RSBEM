function [OK_points] = filter_inside_or_near_boundaries(points,min_distance,Solutions,Meshes, do_waitbar)

%% Filter out points inside body or tail and points within distance tolerance of moving boundaries

% min_distance = 1.5; % 2.5

OK_points = points;
if do_waitbar
    fh = waitbar(0,'points filtering');
end
for phase_ind = 1:length(Solutions.phase)
    if do_waitbar
        waitbar(phase_ind / length(Solutions.phase) );
    end
    
    Mesh = Meshes{phase_ind};
    
    tail_ind = find(strcmp({Mesh.name},'Tail'));
    body_ind = find(strcmp({Mesh.name},'Body'));
    if ~isempty([tail_ind body_ind])
        [body_tail] = combine_Meshes(Mesh, [body_ind tail_ind]); %%leave out transverse, hair sheets since they're 2D
        [is_inside] = inside_mesh(OK_points, body_tail);
    else
        is_inside = false(size(OK_points,1),1);
    end
    
    %         [except_body] = combine_Meshes(   Mesh, setdiff(1:length(Mesh),[tail_ind body_ind])    );  % why is tail_ind here?
    [except_body] = combine_Meshes(   Mesh, setdiff(1:length(Mesh),[body_ind])    );  %this will work even if body isn't there
    
    [corner_elems, flat_verts] = flatten_mesh(except_body);
    [ distances ] = point2trimesh('Faces',corner_elems,'Vertices',flat_verts,'QueryPoints',OK_points, 'Algorithm','parallel');
    too_close = abs(distances) <= min_distance;
    bad_points = (is_inside | too_close);
    
    OK_points(bad_points,:) = [];
    
    
end

if do_waitbar
    close(fh);
end