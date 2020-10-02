function [VL_centroid] = centroid_integrals(Mesh, input)

variable_type = 'centroid';

if length(Mesh) > 1 %multiple submeshes input, use global indices
    tot_elems = Mesh(end).indices.glob.unq_bounds.elem(2 );
else %only one submesh input, use local indices
    tot_elems = Mesh(1).n_elem;
end

rule = 2;  %1st rule has ridiculously innaccurate (too conservative) error estimates, so never converges to the tolerances
[rule_constants.G, rule_constants.Weights, rule_constants.PTS] =  SMPRMS( 2, rule );  %always 2 dimensions and constant rule # for all integrals (2nd arg is rule #)
coder.varsize('VRTS',[2 3 Inf],[true true true]);
VRTS = [0 1 0; 0 0 1];
integrand_constants = []; %placeholder

[vert_range, elem_range] = bounds_from_global_inds(Mesh);
VL_centroid = NaN(tot_elems,18); 

  [i_mesh_elems, local_elems] = global2local(1:tot_elems, Mesh, elem_range,'elem');


parfor (elem_i = 1:tot_elems, input.performance.nthreads)
%    for elem_i = 1:tot_elems

i_mesh_elem = i_mesh_elems(elem_i);
local_elem = local_elems(elem_i);

%     if length(Mesh) > 1 %multiple submeshes
%         [i_mesh_elem, local_elem] = global2local(elem_i, Mesh, elem_range,'elem');
%     else
%         i_mesh_elem = 1;
%         local_elem = elem_i; %we're already looping over local elems in case of one submesh
%     end
    %for elem_i = 1:Mesh.n_elem
    subverts = Mesh(i_mesh_elem).verts(Mesh(i_mesh_elem).elems(local_elem,:),:);
    
    shape_parameters = Mesh(i_mesh_elem).elem_params(:,local_elem);  %[alpha beta gamma]
    
     %scale abs integration tolerance by element area, since big elements
        %should be allowed more total error than small elements
        abstol =  input.accuracy.integration_tol.(variable_type).abstol * Mesh(i_mesh_elem).area(local_elem);
        %don't scale reltol, since reltol should already account for area
        %since the integral itself should scale with area

    [VL, ~, ~, ~] = adsimp( 2, VRTS, 18,  input.accuracy.integration_tol.(variable_type).maxevals, abstol, input.accuracy.integration_tol.(variable_type).reltol ,rule  , rule_constants,subverts,shape_parameters,integrand_constants,'centroid');  
  %VL is a 18 element vector that is the integral of [1/2 x^2  0  0] dot n * hs * phi , [0  1/2 y^2  0] dot n * hs * phi , [0  0  1/2 z^2] dot n * hs * phi   for this element
    VL_centroid(elem_i,:) = VL;

end
