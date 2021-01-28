

% compare velocities on mesh nodes based on U and Omega vs field velocity calculation


% run main_test2 up to first matrix assembly, then solve with sol = A_temp \ RHS
% then calc field velocities at mesh nodes

clear pts
bads = [];
n_pts = 10;
vals = linspace(0,1,n_pts);
c = 0;
for i = 1:Mesh(1).n_elements
    for j1 = 1:n_pts
        xi = vals(j1);
        for j2 = 1:n_pts
            
            eta = vals(j2);
            if eta > (1 - xi) % not filling in an entire square, only a triangle bounded above by eta = 1 - xi
                break
            end
            c = c + 1;
            %             [xi eta]
            nodes = Mesh(1).nodes(Mesh(1).elements(i,:),:);
            pts(c,:) = T6interp(nodes,xi,eta,Mesh(1).shape_parameters(i,:));
%             dists = sqrt(sum((pts(c,:) - nodes).^2,2));
%             if ~any(dists < 1E-10)
%                 bads(end+1) = i;
%             end
        end
    end
end

r = sqrt(sum(pts.^2,2));
[min(r) mean(r) max(r)]
%%
 figure(349); [s,e] = plot_mesh(Mesh(1),2);  e.EdgeAlpha = 0.7;
 hold on; plot_mesh(Mesh(2),1);
 
 
temp_input = assembly_input; %temp_input.accuracy.mesh.integration_tol.stokeslet.abstol = 0.002;
[u_field] = field_velocity(Mesh,Network, s.Vertices, solution, matrix_props,index_mapping,mesh_node_parameters,temp_input);


%%
% [X,Y,Z] = sphere(100);
% sph = surface(X*Metadata.geom.sphererad,Y*Metadata.geom.sphererad,Z*Metadata.geom.sphererad);
% sph.FaceAlpha = 0.4; sph.EdgeAlpha = 0.1;



der_swimmer = solution(end-6:end);

clear u_rigid
for i = 1:size(s.Vertices,1)
    
    u_rigid(i,:) = der_swimmer(1:3) + crossprod(der_swimmer(4:6), s.Vertices(i,:)' - Mesh(2).refpoints);
    
end
%%

for i = 1:size(Mesh(2).nodes,1)
    ii = size(Mesh(1).nodes,1) + i;
    
    u_rigid(ii,:) = der_swimmer(1:3) + crossprod(  der_swimmer(4:6) - der_swimmer(7)*Mesh(2).orientation(:,1)     , Mesh(2).nodes(i,:)' - Mesh(2).refpoints);
    
end

%%
speed_rigid = sqrt(sum(u_rigid.^2,2));
 speed_field = sqrt(sum(u_field.^2,2));
 diffs = (speed_rigid - speed_field);
 reldiffs = diffs / mean(speed_rigid);
 
 [min(reldiffs)  mean(reldiffs)  max(reldiffs)]
 

 s.FaceVertexCData = reldiffs;
s.FaceColor = 'interp';
colorbar;

%%
% ans =
%     1.557e-14   0.00062909     0.046998    % nodes about to enter
% ans =
%    3.9482e-13   0.00042563    0.0034241   % t = 0
% ans =
%     3.786e-13   0.00042264    0.0039364   % nodes about to enter, but setting E = 0 so network doesn't matter here
%%

figure;  plot_mesh(Mesh,2); light;


hold on;  try, delete(fq5); end; factor = 1; fq5 = quiver3([Mesh(1).nodes(:,1); Mesh(2).nodes(:,1)],[Mesh(1).nodes(:,2); Mesh(2).nodes(:,2)],[Mesh(1).nodes(:,3); Mesh(2).nodes(:,3)],...
    [u_rigid(:,1)],[u_rigid(:,2); ],[u_rigid(:,3);]);  fq5.Color = 'k'; fq5.LineWidth = 2; fq5.MaxHeadSize = 10;

hold on;  try, delete(fq4); end; factor = 1; fq4 = quiver3([Mesh(1).nodes(:,1); Mesh(2).nodes(:,1)],[Mesh(1).nodes(:,2); Mesh(2).nodes(:,2)],[Mesh(1).nodes(:,3); Mesh(2).nodes(:,3)],...
    u_field(:,1),u_field(:,2),u_field(:,3));  fq4.Color = 'r'; fq4.LineWidth = 2; fq4.MaxHeadSize = 10;



%
% hold on;  try, delete(fq5); end; factor = 1; fq5 = quiver3([Mesh(1).nodes(:,1); Mesh(2).nodes(:,1)],[Mesh(1).nodes(:,2); Mesh(2).nodes(:,2)],[Mesh(1).nodes(:,3); Mesh(2).nodes(:,3)],...
%     [u_rigid{1}(:,1); u_rigid{2}(:,1)],[u_rigid{1}(:,2); u_rigid{2}(:,2)],[u_rigid{1}(:,3); u_rigid{2}(:,3)]);  fq5.Color = 'k'; fq5.LineWidth = 2; fq5.MaxHeadSize = 10;
%
% hold on;  try, delete(fq4); end; factor = 1; fq4 = quiver3([Mesh(1).nodes(:,1); Mesh(2).nodes(:,1)],[Mesh(1).nodes(:,2); Mesh(2).nodes(:,2)],[Mesh(1).nodes(:,3); Mesh(2).nodes(:,3)],...
%     u_field(:,1),u_field(:,2),u_field(:,3));  fq4.Color = 'r'; fq4.LineWidth = 2; fq4.MaxHeadSize = 10;



