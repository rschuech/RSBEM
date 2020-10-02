function [signed_flow] = integrate_flow(Mesh, field_vel, input)
% input Mesh should be Mesh_clearance

tot_elems = 0;
for i = 1:length(Mesh)
    tot_elems = tot_elems + Mesh(i).n_elem;
end

variable_type = 'field_vel';

rule = 2;  %1st rule has ridiculously inaccurate (too conservative) error estimates, so never converges to the tolerances
[rule_constants.G, rule_constants.Weights, rule_constants.PTS] =  SMPRMS( 2, rule );  %always 2 dimensions and constant rule # for all integrals (2nd arg is rule #)
coder.varsize('VRTS',[2 3 Inf],[true true true]);
VRTS = [0 1 0; 0 0 1];
% integrand_constants.field_vel = NaN(6,3);

signed_flow = NaN(tot_elems,2);

start = 1;  %start ind for where to slot in current mesh's VL_unity_temp results into VL_unity
for i_mesh = 1:length(Mesh)
    signed_flow_temp = NaN(Mesh(i_mesh).n_elem,2);
    

    parfor (local_elem = 1:Mesh(i_mesh).n_elem, input.performance.nthreads)  %counts over local elem inds, so elements shared over multiple submeshes would be counted multiple times, but unclear what one "should" do in that very weird situation....
%         for local_elem = 1:Mesh(i_mesh).n_elem
        

        subverts = Mesh(i_mesh).verts(Mesh(i_mesh).elems(local_elem,:),:);
        
        shape_parameters = Mesh(i_mesh).elem_params(:,local_elem);  %[alpha beta gamma]
        %scale abs integration tolerance by element area, since big elements
        %should be allowed more total error than small elements
        abstol =  input.accuracy.integration_tol.(variable_type).abstol * Mesh(i_mesh).area(local_elem);
        %don't scale reltol, since reltol should already account for area
        %since the integral itself should scale with area
        integrand_constants = [];  
        integrand_constants.field_vel = field_vel(Mesh(i_mesh).elems(local_elem,:),:); % don't have field_vel stored for multiple submeshes currently
        [VL, ~, ~, ~] = adsimp( 2, VRTS, 2,  input.accuracy.integration_tol.(variable_type).maxevals, abstol, input.accuracy.integration_tol.(variable_type).reltol ,rule  , rule_constants,subverts,shape_parameters,integrand_constants,'signed_flow');
       
        signed_flow_temp(local_elem,:) = VL; % first is positive flow, 2nd is negative flow
        
    end
    
    signed_flow(start:(start + Mesh(i_mesh).n_elem - 1),:) = signed_flow_temp;
    start = start + Mesh(i_mesh).n_elem;
    
end


