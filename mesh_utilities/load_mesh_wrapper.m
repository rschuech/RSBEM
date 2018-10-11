function [Mesh, Metadata, rand_inds] = load_mesh_wrapper(meshname,input,rand_inds,diffs)


if isempty(diffs)  %only differences were above fields, so we must have just done forced and now freeswim, or just done freeswim and now forced, everything else the same
    [Mesh, Metadata] = load_mesh(meshname, input, rand_inds); %use rand_inds from before so that randomized Mesh struct is identical to before
    
else %something else changed, so cannot reuse any part of A matrix (or this is the first sweep iteration)
    [Mesh, Metadata] = load_mesh(meshname, input); %by convention, index 1 used for body mesh and geometry
    
end

if isfield(Metadata.mesh,'rand_inds') %should be there, if we are randomizing verts
    rand_inds = Metadata.mesh.rand_inds;  %store for next sweep iter
else %use placeholder since it won't actually be used  (this only happens when input.performance.randomize_verts == false)
    rand_inds = [];
end