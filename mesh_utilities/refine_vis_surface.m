function [Surface_vis] = refine_vis_surface(Mesh, n_refines, surfun)

% levels of refinement: n_refines = 1 uses existing midpoint nodes to split each triangle into 4;
% n_refines > 1 further splits each of those triangles into 4 new triangles,
% (n_refines - 1) more times

%n_refines can either be a scalar, in which case it is used for all
%submeshes, or a vector equal to the number of submeshes = length(Mesh)

%returns Surface_vis, a refined mesh comprised of flat 3 vertex triangles that
%can be plotted with patch()

% surfun is an optional cell array (# cells = length(Mesh), each cell = n_vert x N components) of values defined over the surface such as
% traction to be refined

if isscalar(n_refines)
    n_refines = repmat(n_refines,1,length(Mesh)); %expand a scalar input of n_refines
end


Surface_vis = Mesh;  %initialize as input, all fields will be copied
% [Surface_vis.locals] = deal(NaN);   %initialize new fields so they exist
% [Surface_vis.xi_eta] = deal(NaN);
% [Surface_vis.xi_eta_elem] = deal(NaN);

for i = 1:length(Surface_vis)
    
    % Surface_vis.verts = Mesh.verts;
    Surface_vis(i).elems = NaN(4*Mesh(i).n_elements,3);
    Surface_vis(i).locals = NaN(4*Mesh(i).n_elements,6);
    Surface_vis(i).xi_eta = NaN(Mesh(i).n_nodes,2);
    Surface_vis(i).xi_eta_elem = NaN(Mesh(i).n_nodes,1);
    % Surface_vis.n_vert = NaN;
    % Surface_vis.n_elem = NaN;
    if nargin == 3
        Surface_vis(i).surfun = surfun{i};
    end
    
end


[Surface_vis] = refine_vis_surface_core(Mesh, n_refines, Surface_vis);


% recompute global indices
%Surface_vis = global_inds(Surface_vis);  this is broken after changing
%mesh conventions for dino, would need to fix if ever needed...


 