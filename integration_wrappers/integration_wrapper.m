function [VL_out] = integration_wrapper(Mesh, integral_type, abstol_unscaled, reltol, maxevals, nthreads, triangle_integration, varargin)
%used to integrate number 1 over the surface
%outputs VL_out, (n_elem x 6) containing integrals of just (1) * hS * phi
% n_elem used to be total # elements I guess

% updated to output a cell array, cells correspond to each submesh, each cell is (Mesh(i).n_elements x 6) containing integrals of just (1) * hS * phi

% if length(Mesh) > 1 %multiple submeshes input, use global indices
%     tot_elems = Mesh(end).indices.glob.unq_bounds.elem(2 );
% else %only one submesh input, use local indices
%     tot_elems = Mesh(1).n_elem;
% end

% tot_elems = 0;
% for i = 1:length(Mesh)
%     tot_elems = tot_elems + Mesh(i).n_elem;
% end
if isempty(varargin)
    refpoint = NaN(3,1);
else
    refpoint = varargin{1};
end

% rule = 2;  %1st rule has ridiculously inaccurate (too conservative) error estimates, so never converges to the tolerances
% [rule_constants.G, rule_constants.Weights, rule_constants.PTS] =  SMPRMS( 2, rule );  %always 2 dimensions and constant rule # for all integrals (2nd arg is rule #)
% coder.varsize('VRTS',[2 3 Inf],[true true true]);
% VRTS = [0 1 0; 0 0 1];
% integrand_constants = []; %placeholder
%[~, elem_range] = bounds_from_global_inds(Mesh);
% VL_out = NaN(sum([Mesh.n_elements]),6);
VL_out = cell(1,length(Mesh));

switch integral_type
    case {"unity", "flux"}
        n_components = 6;
    case {"centroid", "moments"}
        n_components = 18;
end
% start = 1;  %start ind for where to slot in current mesh's VL_out_temp results into VL_out
for i_mesh = 1:length(Mesh)
    
    VL_out_temp = NaN(Mesh(i_mesh).n_elements,n_components);
%     [i_meshs, local_elems] = global2local(1:tot_elems, Mesh, elem_range,'elem');
    
    parfor (local_elem = 1:Mesh(i_mesh).n_elements, nthreads)  %counts over local elem inds, so elements shared over multiple submeshes would be counted multiple times, but unclear what one "should" do in that very weird situation....
%         for local_elem = 1:Mesh(i_mesh).n_elements
            %for elem_i = 1:tot_elems
        
        %     i_mesh = i_meshs(elem_i);
        %  local_elem = local_elems(elem_i);
        
        %     if length(Mesh) > 1 %multiple submeshes
        %         [i_mesh, local_elem] = global2local(elem_i, Mesh, elem_range,'elem');
        %     else
        %         i_mesh = 1;
        %         local_elem = elem_i; %we're already looping over local elems in case of one submesh
        %     end
         integrand_constants = struct('eps2',NaN,...
        'x_col',NaN(3,1),...
        'element_nodes',Mesh(i_mesh).nodes(Mesh(i_mesh).elements(local_elem,:),:),...
        'shape_parameters',Mesh(i_mesh).shape_parameters(local_elem,:),...
        'integral_type',integral_type,...
        'DL_singularity_removal',false,... % doesn't matter for single layer integral, will overwrite later for double layer integral
        'coll_node_ind', NaN,...
        'refpoint', refpoint);
        
%         integrand_constants.element_nodes = Mesh(i_mesh).nodes(Mesh(i_mesh).elements(local_elem,:),:);
        
%         integrand_constants.shape_parameters = Mesh(i_mesh).shape_parameters(:,local_elem);  %[alpha beta gamma]
        %scale abs integration tolerance by element area, since big elements
        %should be allowed more total error than small elements
         abstol_scaled =  abstol_unscaled * Mesh(i_mesh).area(local_elem);
        %don't scale reltol, since reltol should already account for area
        %since the integral itself should scale with area
        
        [VL, ~, ~, ~] = adsimp( 2, triangle_integration.reference_nodes, n_components, maxevals, abstol_scaled, reltol  , triangle_integration,integrand_constants);
        %VL is a 6 element vector that is the integral of hs * phi for this element
        VL_out_temp(local_elem,:) = VL;
        
    end
    
%     VL_out(start:(start + Mesh(i_mesh).n_elements - 1),:) = VL_out_temp;
%     start = start + Mesh(i_mesh).n_elem;

VL_out{i_mesh} = VL_out_temp;
    
end

%to integrate force, really just need to integrate 1 (i.e. the number 1) over the surface,
%then finally matrix multiply by f later on
