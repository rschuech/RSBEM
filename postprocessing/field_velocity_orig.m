function [field_vel] = field_velocity(points, f, Mesh, field_input)
% points = N x 3
% field_vel = N x 3


rule = 2;  %1st rule has ridiculously innaccurate (too conservative) error estimates, so never converges to the tolerances
[rule_constants.G, rule_constants.Weights, rule_constants.PTS] =  SMPRMS( 2, rule );  %always 2 dimensions and constant rule # for all integrals (2nd arg is rule #)
coder.varsize('VRTS',[2 3 Inf],[true true true]);
VRTS = [0 1 0; 0 0 1];
%%

fmat = reshape(f,3,[]);  %row is x y z traction, col is global vert

% U = NaN(size(X));  V = U;  W = U;
field_vel = NaN(size(points));

if length(Mesh) > 1 %multiple submeshes input, use global indices
    tot_elems = Mesh(end).indices.glob.unq_bounds.elem(2 );
else %only one submesh input, use local indices
    tot_elems = Mesh(1).n_elem;
end

[vert_range, elem_range] = bounds_from_global_inds(Mesh);

parfor field_i = 1:size(points,1)
    
    %do boundary integral over all mesh elements to get field velocity at this field point
    VL_fieldvel = NaN(tot_elems,3);
    
    for elem_i = 1:tot_elems  %elem_i is a global index
        %                  elem_i / tot_elems
        if length(Mesh) > 1 %multiple submeshes
            [i_mesh_elem, local_elem] = global2local(elem_i, Mesh, elem_range, 'elem');
        else
            i_mesh_elem = 1;
            local_elem = elem_i; %we're already looping over local elems in case of one submesh
        end
        
        vertinds = Mesh(i_mesh_elem).elems(local_elem,:);
        subverts = Mesh(i_mesh_elem).verts(vertinds,:);
        global_vertinds = Mesh(i_mesh_elem).indices.glob.vert(vertinds);
        
        sub_f = fmat(:,global_vertinds)';
        shape_parameters = Mesh(i_mesh_elem).elem_params(:,local_elem);  %[alpha beta gamma]
        
        abstol =  field_input.accuracy.integration_tol.field_vel.abstol * Mesh(i_mesh_elem).area(local_elem);
        %don't scale reltol, since reltol should already account for area
        %since the integral itself should scale with area
        
        %the important part:  adaptively compute boundary integrals of reg stokeslet function over
        %current element, with respect to current collocation point
        %integrand_constants = struct('eps2',assembly_input.accuracy.eps2,'x_col',Mesh(i_mesh_vert).verts(local_vert,:)');
        integrand_constants = struct('f',sub_f,'x_field',points(field_i,:)', 'eps2',field_input.accuracy.eps2);
        %         integrand_constants = [];
        %         integrand_constants.f = sub_f;  %like subverts, except the x y z traction components instead of coordinates at each vertex of the element
        %         integrand_constants.x_field = [X(field_i) Y(field_i) Z(field_i)]';
        %         integrand_constants.eps2 = input.accuracy.eps2;
        %     if assembly_input.performance.debug_mode
        [VL, AE, NV, FL] = adsimp( 2, VRTS, 3,  field_input.accuracy.integration_tol.field_vel.maxevals,abstol, field_input.accuracy.integration_tol.field_vel.reltol ,rule  , rule_constants ,subverts,shape_parameters,integrand_constants,'S_dot_f');
        VL_fieldvel(elem_i,:) = VL;
        %         %
    end
    
    field_vel(field_i,:) = sum(VL_fieldvel,1);
    %U(field_i) = temp(1);  V(field_i) = temp(2);  W(field_i) = temp(3);
    
end

field_vel = field_vel * field_input.constants.multfactor / field_input.constants.mu; %scale by term in front of boundary integral equation