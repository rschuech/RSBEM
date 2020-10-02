field_input.accuracy.integration_tol.field_vel.abstol = 1E-5;  % apparently too small for body-centered-tail but 1E-5, 1E-6 give nearly same results for that one so prolly fine
% field_input.accuracy.integration_tol.field_vel.abstol = 1E-6;
%field_input.accuracy.integration_tol.field_vel.abstol = 1E-5;

field_input.accuracy.integration_tol.field_vel.reltol = 0;
field_input.accuracy.integration_tol.field_vel.maxevals = Inf;
field_input.accuracy.eps2 = input.accuracy.eps2;
field_input.constants.multfactor = input.constants.multfactor;
field_input.constants.mu = input.constants.mu;

field_input.performance.nthreads = 20;



%% 


[clearance_mesh, clearance_metadata] = load_mesh('C:\Users\rudi\Desktop\RD\clearance rate\ClearanceRateMesh.dat',[],[],'mesh');

factor = 1.5;

Mesh_c = clearance_mesh;


[Mesh_c] = global_inds(Mesh_c);
[Mesh_c] = renumber_Mesh(Mesh_c);

clear temp_input
temp_input.performance.nthreads = 8;
temp_input.accuracy.integration_tol.area.abstol = 1E-6;
temp_input.accuracy.integration_tol.area.reltol = 1000;
temp_input.accuracy.integration_tol.area.maxevals = Inf;

temp_input.accuracy.integration_tol.centroid.abstol = 1E-6;
temp_input.accuracy.integration_tol.centroid.reltol = 1000;
temp_input.accuracy.integration_tol.centroid.maxevals = Inf;

temp_input.accuracy.integration_tol.volume.abstol = 1E-6;
temp_input.accuracy.integration_tol.volume.reltol = 1000;
temp_input.accuracy.integration_tol.volume.maxevals = Inf;

[Mesh_c] = store_mesh_constants_wrapper(Mesh_c, temp_input);

shift = Mesh_c.Centroid;
Mesh_c = shiftMesh(Mesh_c,-shift);
Mesh_c.verts = Mesh_c.verts * factor;  Mesh_c.area = Mesh_c.area*factor^2;  Mesh_c.centroid = Mesh_c.centroid * factor;  Mesh_c.Volume = Mesh_c.Volume*factor^3;
Mesh_c = shiftMesh(Mesh_c,shift);


%%
clear mL_per_hr bodyvols_per_hr
for phase_ind = 1:length(Solutions.f)
    phase_ind
% phase_ind = 16 ;
%         points = Mesh_c.centroid'; % faster than using all verts, just to test if body frame speeds are zero or close enough
points = Mesh_c.verts ; % need to use verts to integrate_flow

[Mesh, Metadata, matrix_props, Metadata_rand_inds] = load_dino_mesh(Solutions.phase(phase_ind), Mesh_files, input, Solutions.rand_inds{phase_ind}); 
if isfield(field_input,'performance'), temp = rmfield(field_input,'performance'); else, temp = field_input; end
[field_vel] = field_velocity_mexed(points, Solutions.f{phase_ind}, Mesh, temp);  % Mesh is the full mesh correspondng to solution f

speed = sqrt(sum(field_vel.^2,2));

omega_cross_r = cross(repmat(Solutions.Omega(phase_ind,:),size(points,1),1)  ,   points - repmat(Mesh(1).refpoints(:,1)',size(points,1),1)   );
field_vel_bodyframe = field_vel - repmat(Solutions.U(phase_ind,:),size(field_vel,1),1) - omega_cross_r;

speed_bodyframe = sqrt(sum(field_vel_bodyframe.^2,2));

field_input.performance.nthreads = 20;
signed_flow = integrate_flow(Mesh_c, field_vel_bodyframe, field_input);  % flow in and flow out of Mesh_c

%  26615     -26587

signed_Flow = sum(signed_flow);

sum(signed_Flow) / min(abs(signed_Flow)) * 100;

mL_per_hr(phase_ind,:) = signed_Flow * 3.6E-9;  % micron^3 / s to mL / hr as in Lasse pape

bodyvols_per_hr(phase_ind,:) = signed_Flow / Mesh(1).Volume * 3600; % 10^5 body volumes / hr expected for flagellates according to Kiorboe

end

% nanmean(mL_per_hr)
% nanmean(bodyvols_per_hr)

  interpolant = spline([Solutions.phase ]',mL_per_hr');
  1/diff(phase_bounds) * integral(@(x) fnval(interpolant, x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12,'ArrayValued',true)
  
   interpolant = spline([Solutions.phase ]',bodyvols_per_hr');
  1/diff(phase_bounds) * integral(@(x) fnval(interpolant, x), phase_bounds(1),phase_bounds(2),'reltol',1E-9,'abstol',1E-12,'ArrayValued',true)

% body transverse tail 2 Coplanar 1.5 Normal (swims at 157 um/s)
% 6.7E-4 mL / hr
% 3.15E4 body vol / hr

% body centered tail (swims at 21.9 um/s)
% 9.5E-5 mL / hr (scaled up by 157/22 = 6.8E-4)
% 4.5E3 body vol / hr (scaled up by 157/22 = 3.2E4)
%%
figure(344); clf; [s,e] = plot_mesh(Mesh);  [sc,ec] = plot_mesh(Mesh_c,0);  sc.FaceAlpha = 0.4;  hold on;
set(s,'facecolor',repmat(0.25,1,3));
% qc = quiver3(Mesh_c.verts(:,1),Mesh_c.verts(:,2),Mesh_c.verts(:,3),field_vel_bodyframe(:,1),field_vel_bodyframe(:,2),field_vel_bodyframe(:,3),4,'k','linewidth',2);
flow = sum( signed_flow, 2 );  filter = flow < 0;
qc = quiver3(Mesh_c.centroid(1,filter)',Mesh_c.centroid(2,filter)',Mesh_c.centroid(3,filter)',field_vel_bodyframe(filter,1),field_vel_bodyframe(filter,2),field_vel_bodyframe(filter,3),1.5,'k','linewidth',1.5);

sc.FaceVertexCData = sum( signed_flow, 2 );
sc.FaceColor = 'flat';
view([140 45]);
axis off
%% for testing with forced (bacteria) runs
% [field_vel] = field_velocity_mexed(points, solutions.rx.f, Mesh, rmfield(field_input,'performance'));  % Mesh is the full mesh correspondng to solution f
% U = [0 0 0];  Omega = [1 0 0];
% %  r_cross_omega = cross( points - repmat(Mesh(2).refpoints(:,1)',size(points,1),1)   ,  repmat(Omega,size(points,1),1)  );
%   omega_cross_r = cross(  repmat(Omega,size(points,1),1) , points - repmat(Mesh(2).refpoints(:,1)',size(points,1),1)     );
% field_vel_bodyframe = field_vel - repmat(U,size(field_vel,1),1) - omega_cross_r;
% speed_bodyframe = sqrt(sum(field_vel_bodyframe.^2,2));
% s.CData = speed_bodyframe;

export_fig(gcf,'C:\Users\rudi\Desktop\RD\pape main results\figures\clearance rate all vectors.png','-r300','-transparent')

%%
phase_ind = 17;
[Mesh, Metadata, matrix_props, Metadata_rand_inds] = load_dino_mesh(Solutions.phase(phase_ind), Mesh_files, input, Solutions.rand_inds{phase_ind}); 

n_unq_verts = 0;
for i = 1:length(Mesh)
    n_unq_verts = max(n_unq_verts, Mesh(i).indices.glob.unq_bounds.vert(2));
end
n_unq_verts*3

size( Solutions.f{phase_ind} ,1)
