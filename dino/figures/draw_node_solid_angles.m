i = 1;

 angles = Mesh(i).solid_angle;  fractions = angles / 2/pi;
 
 
figure(234)
clf;  [s,e] = plot_mesh(Mesh(i),2);  e.EdgeAlpha = 0.4; hold on

try, delete(sc), end;  filter = fractions > 1.05;  sc = scatter3(Mesh(i).verts(filter,1),Mesh(i).verts(filter,2),Mesh(i).verts(filter,3),30,'k','filled');
view([0 0])
axis tight
axis off