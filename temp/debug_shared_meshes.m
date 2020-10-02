find ( Mesh(2).indices.glob.vert < Mesh(2).indices.glob.unq_bounds.vert(1) | Mesh(2).indices.glob.vert > Mesh(2).indices.glob.unq_bounds.vert(2))

shat = ismember(Mesh(2).elems, inds);
% verts that are on transverse and shared with body

shared_glob = Mesh(2).indices.glob.vert( Mesh(2).indices.glob.vert < Mesh(2).indices.glob.unq_bounds.vert(1) | Mesh(2).indices.glob.vert > Mesh(2).indices.glob.unq_bounds.vert(2));



% verts that are on body that are shared with transverse
body_vert_inds = find( ismember(Mesh(1).indices.glob.vert,shared_glob) );

%elems that are on body and share verts with transverse
body_elem_inds = find( any( ismember(Mesh(1).elems, body_vert_inds) , 2) );


body = Mesh(1);
body.elems = body.elems(body_elem_inds,:);
body.elem_params = body.elem_params(:,body_elem_inds);
body.n_elem = size(body.elems,1);



Surface_vis = refine_vis_surface(body,4);
Edge_vis = refine_vis_edges(body,15);  


for i = 1:length(Surface_vis)
    s(i) = patch('faces',Surface_vis(i).elems,'vertices',Surface_vis(i).verts,'edgecolor','none','facecolor',colors{i},'facealpha',facealpha);
    hold on
    e(i) = patch('faces',Edge_vis(i).elems,'vertices',Edge_vis(i).verts,'facecolor','none','edgecolor','k','edgealpha',1,'linewidth',2);
end

hold off
axis equal; grid on;  box on;
xlabel('x');  ylabel('y');  zlabel('z');