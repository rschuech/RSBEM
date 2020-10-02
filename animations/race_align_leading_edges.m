
% t = 0;
clear leading_x

for d = 1:length(Meshes)
    
%     interpolant = Interpolants{d};
    Mesh = Meshes{d};
    
%     y = ppval(interpolant, t);
    
%     [Mesh(2)] = rotateMesh(Mesh(2), [0 0 y(7)]' );  %rotate tail around x axis
    %         y(3) = y(3) + shifts(d);
    
    % was going to fix this part for pole2pole tails but actually don't see
    % why it is needed at all
    
%     temp = Mesh_temp(2).refpoints(:,1);
%     Mesh_temp(2) = shiftMesh(Mesh_temp(2), -Mesh_temp(2).refpoints(:,1));  % shift so that refpoint is at origin to do rotation around motor axis
%     Mesh_temp(2) = rotateMesh(Mesh_temp(2), Mesh_temp(2).orientation(:,1) );
%     Mesh_temp(2) = shiftMesh(Mesh_temp(2), temp);  % shift back to where it was
%     
    
    
%     Mesh = move_Mesh(Mesh,y);
    
    if align_paths
        Mesh = shiftMesh(Mesh,shift(d,:));
        Mesh = rotateMesh(Mesh,angle(d,:),rotvec(d,:));
    end
    
    Mesh = shiftMesh(Mesh,[0 0 shifts(d)]);
    
    
    leading_x(d) = max(Mesh(1).verts(:,1));
    
end




%         race_shifts = max(leading_x) - leading_x; % aligns leading edges
%         at arbitrary x
race_shifts =  - leading_x; % aligns leading edges at x = 0 in final video



