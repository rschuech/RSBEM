function [VL_flux] = integrate_flux(Mesh, input, variable_type)
%integrates flux of vector field F through surface
%(used via divergence theorem to get mesh volume)

tot_elems = 0;
for i = 1:length(Mesh)
    tot_elems = tot_elems + Mesh(i).n_elem;
end


rule = 2;  %1st rule has ridiculously innaccurate (too conservative) error estimates, so never converges to the tolerances
[rule_constants.G, rule_constants.Weights, rule_constants.PTS] =  SMPRMS( 2, rule );  %always 2 dimensions and constant rule # for all integrals (2nd arg is rule #)
coder.varsize('VRTS',[2 3 Inf],[true true true]);
VRTS = [0 1 0; 0 0 1];
integrand_constants = []; %placeholder
%[~, elem_range] = bounds_from_global_inds(Mesh);
VL_flux = NaN(tot_elems,6);

% [i_meshs, local_elems] = global2local(1:tot_elems, Mesh, elem_range,'elem');
start = 1;  %start ind for where to slot in current mesh's VL_unity_temp results into VL_unity
for i_mesh = 1:length(Mesh)
    VL_flux_temp = NaN(Mesh(i_mesh).n_elem,6);
    
    parfor (local_elem = 1:Mesh(i_mesh).n_elem, input.performance.nthreads)
        %for local_elem = 1:tot_elems
        
        %     i_mesh = i_meshs(local_elem);
        % local_elem = local_elems(local_elem);
        
        %     if length(Mesh) > 1 %multiple submeshes
        %         [i_mesh, local_elem] = global2local(local_elem, Mesh, elem_range,'elem');
        %     else
        %         i_mesh = 1;
        %         local_elem = local_elem; %we're already looping over local elems in case of one submesh
        %     end
        %for local_elem = 1:Mesh.n_elem
        subverts = Mesh(i_mesh).verts(Mesh(i_mesh).elems(local_elem,:),:);
        
        shape_parameters = Mesh(i_mesh).elem_params(:,local_elem);  %[alpha beta gamma]
        
        %scale abs integration tolerance by element area, since big elements
        %should be allowed more total error than small elements
        abstol =  input.accuracy.integration_tol.volume.abstol * Mesh(i_mesh).area(local_elem);
        %don't scale reltol, since reltol should already account for area
        %since the integral itself should scale with area
        
        [VL, ~, ~, ~] = adsimp( 2, VRTS, 6,  input.accuracy.integration_tol.(variable_type).maxevals, abstol, input.accuracy.integration_tol.(variable_type).reltol ,rule  , rule_constants,subverts,shape_parameters,integrand_constants,'flux');
        
        VL_flux_temp(local_elem,:) = VL;
        
    end
    
    VL_flux(start:(start + Mesh(i_mesh).n_elem - 1),:) = VL_flux_temp;
    start = start + Mesh(i_mesh).n_elem;
    
end
