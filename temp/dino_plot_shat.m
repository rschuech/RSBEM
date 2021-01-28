% load a dump with a Mesh here

Mesh_p = Mesh([2 3 4  6   ]);

% Mesh_p(2) = shiftMesh(Mesh_p(2),geom.tail.shift);
% 
figure(838); clf;  [s,e] = plot_mesh(Mesh_p,[2 2 2 2 2 2 2 2]);  view([-20 25]);
view([-20 30]);
view([-90 0]);
% li = light('position',[0 -1 1]);
li = light('position',[-1 0 0]); li2 = light('position',[0 0 1]);
xlim([-11.5 48]); ylim([-24 24]); zlim([-24 24]); % most cases
% xlim([-32 48]); ylim([-24 24]); zlim([-24 24]);  % centered tail
axis off

% export_fig(gcf,'C:\Users\rudi\Desktop\RD\pape main results\figures\Body Centered-Tail.png','-r300','-transparent','-nocrop');