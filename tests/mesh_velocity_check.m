

% compare velocities on mesh nodes based on U and Omega vs field velocity calculation


% run main_test2 up to first matrix assembly, then solve with sol = A_temp \ RHS
% then calc field velocities at mesh nodes

der_swimmer = sol(end-6:end);

clear u_rigid
for i = 1:size(Mesh(1).nodes,1)
    
    u_rigid(i,:) = der_swimmer(1:3) + crossprod(der_swimmer(4:6), Mesh(1).nodes(i,:)' - Mesh(2).refpoints);
    
end


for i = 1:size(Mesh(2).nodes,1)
    ii = size(Mesh(1).nodes,1) + i;
    
    u_rigid(ii,:) = der_swimmer(1:3) + crossprod(  der_swimmer(4:6) - der_swimmer(7)*Mesh(2).orientation(:,1)     , Mesh(2).nodes(i,:)' - Mesh(2).refpoints);
    
end

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



