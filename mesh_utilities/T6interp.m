function [x_interp, phi, hS,  n, e] = T6interp(nodes,xi,eta,shape_parameters)
%nodes = [x1 y1 z1; x2 y2 z2; ... x6 y6 z6] for 6 nodes of a quadratic
%triangle  (6 x 3)

%Pozrikidis page 122

%following definitions may be sped up with simplified formulas from Pozrikidis
%textbook if middle edge nodes are always in the exact middle
% alpha = 1/(1+norm(nodes(4,:)-nodes(2,:))/norm(nodes(4,:)-nodes(1,:)))
% beta = 1/(1+norm(nodes(6,:)-nodes(3,:))/norm(nodes(6,:)-nodes(1,:)))
% gamma = 1/(1+norm(nodes(5,:)-nodes(2,:))/norm(nodes(5,:)-nodes(3,:)))

nodes = nodes';  % now 3 x 6

alpha = shape_parameters(1); % scalar
beta = shape_parameters(2);
gamma = shape_parameters(3);

phi = NaN(6,1);  


phi(2) = 1/(1-alpha)*xi*(xi-alpha+(alpha-gamma)/(1-gamma)*eta);
phi(3) = 1/(1-beta)*eta*(eta-beta+(beta+gamma-1)/gamma*xi);
phi(4) = 1/alpha/(1-alpha)*xi*(1-xi-eta);
phi(5) = 1/gamma/(1-gamma)*xi*eta;
phi(6) = 1/beta/(1-beta)*eta*(1-xi-eta);
phi(1) = 1-phi(2)-phi(3)-phi(4)-phi(5)-phi(6);

x_interp = nodes * phi;  % (3 x 6)  *  (6 x 1) = (3 x 1)

if nargout > 2
    
    %derivatives
    phi_xi = NaN(6,1);  phi_eta = NaN(6,1);  % (1 x 6)
    
    
    phi_xi(2) = 1/(1-alpha) * ( 2*xi - alpha + (alpha-gamma)/(1-gamma)*eta);
    phi_xi(3) = 1/(1-beta)*eta*(beta+gamma-1)/gamma;
    phi_xi(4) = 1/alpha/(1-alpha)*(1-2*xi-eta);
    phi_xi(5) = 1/gamma/(1-gamma)*eta;
    phi_xi(6) = -1/beta/(1-beta)*eta;
    phi_xi(1) = - phi_xi(2) - phi_xi(3) - phi_xi(4) - phi_xi(5) - phi_xi(6);
    
    
    phi_eta(2) = 1/(1-alpha)*xi*(alpha-gamma)/(1-gamma);
    phi_eta(3) = 1/(1-beta)*(2*eta - beta + (beta+gamma-1)/gamma*xi);
    phi_eta(4) = -1/alpha/(1-alpha)*xi;
    phi_eta(5) = 1/gamma/(1-gamma)*xi;
    phi_eta(6) = 1/beta/(1-beta)*(1-xi-2*eta);
    phi_eta(1) = - phi_eta(2) - phi_eta(3) - phi_eta(4) - phi_eta(5) - phi_eta(6);
    
    
    e_xi = nodes * phi_xi; % (3 x 6) * (6 x 1) = (3 x 1)
    e_eta = nodes * phi_eta;
    
%     crossvec = crossprod(e_xi,e_eta); % (3 x 1)
    crossvec = [e_xi(2)*e_eta(3) - e_xi(3)*e_eta(2); e_xi(3)*e_eta(1) - e_xi(1)*e_eta(3); e_xi(1)*e_eta(2) - e_xi(2)*e_eta(1)]; % (3 x 1)
    
    hS = sqrt(sum(crossvec.^2)); % scalar
    
    
    if nargout >= 4 %also output surface normal
        
        n = crossvec / hS; % (3 x 1) supposed to be unit length
%         n = n./ sqrt(sum(n.^2)); % make sure to be safe....
    end
    
    if nargout >= 5 %also output surface tangents
        
       e = [e_xi e_eta]; % (3 x 1, 3 x 1 = 3 x 2)  not positive these are unit length
       e = e ./ sqrt(sum(e.^2,1));
    end
    
    
    
end