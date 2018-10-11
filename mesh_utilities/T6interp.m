function [x_interp, phi, hS,  n] = T6interp(verts,xi,eta,shape_parameters)
%verts = [x1 y1 z1; x2 y2 z2; ... x6 y6 z6] for 6 verts of a quadratic
%triangle

%Pozrikidis page 122

%following definitions may be sped up with simplified formulas from Pozrikidis
%textbook if middle edge nodes are always in the exact middle
% alpha = 1/(1+norm(verts(4,:)-verts(2,:))/norm(verts(4,:)-verts(1,:)))
% beta = 1/(1+norm(verts(6,:)-verts(3,:))/norm(verts(6,:)-verts(1,:)))
% gamma = 1/(1+norm(verts(5,:)-verts(2,:))/norm(verts(5,:)-verts(3,:)))


alpha = shape_parameters(1);
beta = shape_parameters(2);
gamma = shape_parameters(3);

 phi = NaN(1,6);


phi(2) = 1/(1-alpha)*xi*(xi-alpha+(alpha-gamma)/(1-gamma)*eta);
phi(3) = 1/(1-beta)*eta*(eta-beta+(beta+gamma-1)/gamma*xi);
phi(4) = 1/alpha/(1-alpha)*xi*(1-xi-eta);
phi(5) = 1/gamma/(1-gamma)*xi*eta;
phi(6) = 1/beta/(1-beta)*eta*(1-xi-eta);
phi(1) = 1-phi(2)-phi(3)-phi(4)-phi(5)-phi(6);

x_interp = verts' * phi';

if nargout > 2
    
    %derivatives
   phi_xi = NaN(1,6);  phi_eta = NaN(1,6);
  
   
    phi_xi(2) = 1/(1-alpha) * ( 2*xi - alpha + (alpha-gamma)/(1-gamma)*eta);
    phi_xi(3) = 1/(1-beta)*eta*(beta+gamma-1)/gamma;
    phi_xi(4) = 1/alpha/(1-alpha)*(1-2*xi-eta);
    phi_xi(5) = 1/gamma/(1-gamma)*eta;
    phi_xi(6) = -1/beta/(1-beta)*eta;
    phi_xi(1) = - phi_xi(2) - phi_xi(3) - phi_xi(4) - phi_xi(5) - phi_xi(6);
    
  
    phi_eta(2) = 1/(1-alpha)*xi*(alpha-gamma)/(1-gamma);
    phi_eta(3) = 1/(1-beta)*(2*eta - beta + (beta+gamma-1)/gamma*eta);
    phi_eta(4) = -1/alpha/(1-alpha)*xi;
    phi_eta(5) = 1/gamma/(1-gamma)*xi;
    phi_eta(6) = 1/beta/(1-beta)*(1-xi-2*eta);
    phi_eta(1) = - phi_eta(2) - phi_eta(3) - phi_eta(4) - phi_eta(5) - phi_eta(6);
    
    
    e_xi = verts' * phi_xi';
    e_eta = verts' * phi_eta';
    
    crossvec = crossprod(e_xi,e_eta);
    
    hS = sqrt(sum(crossvec.^2));

        
        if nargout == 4 %also output surface normal
            
            n = crossvec / hS;
        end

    
    
    
    
end