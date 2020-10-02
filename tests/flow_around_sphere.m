clear field_input

field_input.accuracy.integration_tol.field_vel.abstol = 1E-2; %5E-9 % apparently too small for body-centered-tail but 1E-5, 1E-6 give nearly same results for that one so prolly fine
% field_input.accuracy.integration_tol.field_vel.abstol = 1E-6;
%field_input.accuracy.integration_tol.field_vel.abstol = 1E-5;

field_input.accuracy.integration_tol.field_vel.reltol = 0;
field_input.accuracy.integration_tol.field_vel.maxevals = Inf;
field_input.accuracy.eps2 = input.accuracy.eps2;
field_input.constants.multfactor = input.constants.multfactor;
field_input.constants.mu = input.constants.mu;

% points = Mesh.verts;  % points on a larger spherical surface around sphere
points = Mesh.centroid';
% temp = (1:20)';
% points = [ Metadata.geom.radius 0 0; 0.7 0 0; [temp zeros(length(temp),2)]  ];

[field_vel_NS] = field_velocity_mexed(points, solutions.x.f, rmfield(Mesh,{'normals','tangents'}), field_input);

vel = field_vel_NS - [repmat(50, size(points,1),1) zeros(size(points,1),2)];
% [field_vel_NS] = field_velocity_mexed(points, solutions.x.f, (Mesh), field_input);
% U = 50 on surface, U = 0 far     fixed frame
% U = 0 on surface, U = -50 far    moving frame

%%


points = Mesh.verts;  % points on a larger spherical surface around sphere
points = [Mesh.verts(:,3)  Mesh.verts(:,1)  Mesh.verts(:,2)];

f = reshape(solutions.x.f,3,[])'; 


clear p P Theta N
r = sqrt(sum(points.^2,2));
theta = acos( dot( points ./r , repmat([1 0 0],size(points,1),1) , 2 ) );
phi = atan2(points(:,3),points(:,2));

% unit normal vector
N = points ./ r;
%unit tangent vector in theta dir
Theta(:,1) = cos(theta + pi/2);
% t2 = -(C*r1+  RootOf(   (r3^2+r2^2)*_Z^2   +2*C*r1*_Z*r3  +   r2^2*C^2   +C^2*r1^2 -r2^2)   *r3  )/r2
% t3 = RootOf((r3^2+r2^2)*_Z^2+2*C*r1*_Z*r3+r2^2*C^2+C^2*r1^2-r2^2)
clear root
for i = 1:size(points,1)
root(i,:) = real(roots([points(i,3).^2 + points(i,2).^2 , 2*cos(theta(i) + pi/2).*points(i,1).*points(i,3) , (points(i,2).^2 + points(i,1).^2).*cos(theta(i) + pi/2).^2 - points(i,2).^2 ]));
end
root = root(:,2);  % both roots always seem to be the same, sans tiny complex components
for i = 1:size(points,1)
  if  abs( cos(theta(i) + pi/2).*points(i,1) + root(i) .* points(i,3) ) <= 1E-7
      Theta(i,2) = 0;
  else
Theta(i,2) = -(cos(theta(i) + pi/2).*points(i,1) + root(i) .* points(i,3) )./ points(i,2);
  end
end
Theta(:,3) = root;


U = [repmat(50,size(points,1),1) zeros(size(points,1),2)];
UdotN = dot(U,N,2);

U_inf = -50;
R = Metadata.geom.radius;

u_theta.no_slip = U_inf * sin(theta) .* (R^3 + 3*R*r.^2 - 4*r.^3) ./ (4*r.^3);
u_r.no_slip = U_inf * cos(theta) .* (R^3 - 3*R*r.^2 + 2*r.^3) ./ (2*r.^3);
p.no_slip = -3*R*U_inf*field_input.constants.mu * cos( theta) ./ (2*r.^2);
sigma_rr.no_slip = 2*field_input.constants.mu*U_inf*cos(theta).*(3*R./(2*r.^2) - 3*R^3./(2*r.^4) ); % stress on surface with normal r, i.e. the sphere
sigma_rtheta.no_slip = -3*field_input.constants.mu*U_inf*sin(theta)*R^3./(2*r.^4); % stress on surface with normal r, i.e. the sphere

u_theta.free_slip = U_inf * sin(theta) .* (R - 2*r) ./ (2*r);
u_r.free_slip = U_inf * cos(theta) .* (r - R) ./ r;
p.free_slip = -R * U_inf * field_input.constants.mu * cos( theta) ./ r.^2;
sigma_rr.free_slip = 2*field_input.constants.mu*U_inf*cos(theta)*R./r.^2;
sigma_rtheta.free_slip = zeros(size(theta));
% cartesian vector form of stresses
P.free_slip = -p.free_slip .* Mesh.normals; % negate pressure because a positive pressure points opposite the normal
Sigma_r.free_slip = sigma_rr.free_slip .* Mesh.normals; % positive viscous normal stress is outward, same as normal
Sigma_theta.free_slip = sigma_rtheta.free_slip .* Theta; % positive shear stress on sphere and Theta are in same direction
Total_stress.free_slip = P.free_slip + Sigma_r.free_slip + Sigma_theta.free_slip;

U_theta.free_slip = u_theta.free_slip .* Theta;
U_r.free_slip = u_r.free_slip .* points;  %cartesian form of radial vectors are point coords
U_tot.free_slip = U_theta.free_slip + U_r.free_slip  +  U;  % fixed frame

drag.free_slip = 6*pi*R*input.constants.mu * 50  * 2/3;


u_phi = 0;
%x = r*sin(theta)*cos(phi)
%y = r*sin(theta)*sin(phi)
%z = r*cos(theta)

% u_x = dr/dt*sin(theta)*cos(phi) + r*dtheta/dt*cos(theta)*cos(phi) - r*sin(theta)*dphi/dt*sin(phi) 
clear u
for field = ["no_slip","free_slip"]
   % U = 0 on surface, U = -50 far   moving frame (direct from original eqs)
   % U = 50 on surface, U = 0 far   fixed frame (after subtracting U_inf)
   % note my axes don't match usual definition for coord conversions
u.(field)(:,2) = u_r.(field).*sin(theta).*cos(phi) + u_theta.(field).*cos(theta).*cos(phi) - u_phi.*sin(theta).*sin(phi) ;
u.(field)(:,3) = u_r.(field).*sin(theta).*sin(phi) + u_theta.(field).*cos(theta).*sin(phi) + u_phi.*sin(theta).*cos(phi);
u.(field)(:,1) = u_r.(field).*cos(theta) - u_theta.(field).*sin(theta) - U_inf;

end

% from Cortez et al 3D RSM pape:
U_Inf = -U_inf;
clear u_FC
u_FC(:,1) = 3*R*U_Inf/4*(1./r.^3 - R^2./r.^5).*points(:,1).^2 + R*U_Inf/4./r.*(3 + R^2./r.^2);
u_FC(:,2) = 3*R*U_Inf/4*(1./r.^3 - R^2./r.^5).*points(:,2).*points(:,1);
u_FC(:,3) = 3*R*U_Inf/4*(1./r.^3 - R^2./r.^5).*points(:,3).*points(:,1);
p_FC = 3/2*R*U_Inf*points(:,1)./r.^3;
P.no_slip = -p_FC .* Mesh.normals; % negate pressure because a positive pressure points opposite the normal
Sigma_r.no_slip = sigma_rr.no_slip .* Mesh.normals; % positive viscous normal stress is outward, same dir as normal
Sigma_theta.no_slip = sigma_rtheta.no_slip .* Theta;
Total_stress.no_slip = P.no_slip + Sigma_r.no_slip + Sigma_theta.no_slip;

f_FC = 3*input.constants.mu / 2 / R * [U_Inf 0 0]'; % traction
drag.no_slip = 6*pi*R*input.constants.mu * 50;


% no slip checks

% check normal (pressure here) stress  f is force on water, so negate once to get force on
% surface, negate again because a + pressure points in toward surface
% instead of along outward normals  (note, for no slip, viscous normal
% stress should be zero so is not included here but generally should be)
[dot(f,Mesh.normals,2)   p_FC p.no_slip];
% check shear stress   negate f to get force on surface
[ -dot(f, Theta,2)   sigma_rtheta.no_slip];

% free slip checks
[ -dot(f, Theta,2)   sigma_rtheta.free_slip];

[dot(f,Mesh.normals,2)    p.free_slip];