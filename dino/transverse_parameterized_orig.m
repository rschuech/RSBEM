function [pt,vel] = transverse_parameterized(u,v, time, geom)

%this is orig out of paper, oriented vertically
% r0 =  [(geom.R+geom.d+geom.w+geom.rf)*cos(u); ...
%        (geom.R+geom.d+geom.w+geom.rf)*sin(u); ...
%         geom.b*(u)] ...
%     + ...
%     geom.rf* ...
%     [cos(u).*cos(geom.c*u - geom.omega*time); ...
%      sin(u).*cos(geom.c*u - geom.omega*time); ...
%      sin(geom.c*u - geom.omega*time)];
% 
%    pt = v.*r0 + [(geom.R+geom.d)*(1-v).*cos(u); (geom.R+geom.d).*(1-v).*sin(u); (1-v).*geom.b.*u];
   
   
   % this is rotated -90 around y-axis as in Salome
   
   r0 = [ - ( geom.b*(u) + geom.rf * sin(geom.c*u - geom.omega*time) ) ; ...
          (geom.R+geom.d+geom.w+geom.rf)*sin(u) + geom.rf* sin(u).*cos(geom.c*u - geom.omega*time); ...
          (geom.R+geom.d+geom.w+geom.rf)*cos(u) +  geom.rf*cos(u).*cos(geom.c*u - geom.omega*time);];

   pt = v.*r0 + [ -(1-v).*geom.b.*u;  (geom.R+geom.d).*(1-v).*sin(u);   (geom.R+geom.d)*(1-v).*cos(u);];
   
   pt(1) = pt(1) + geom.shift; %translate in x dir
   
   if nargout > 1
        r0_t = [ - ( - geom.omega*geom.rf * cos(geom.c*u - geom.omega*time) ) ; ...
            geom.omega*geom.rf* sin(u).*sin(geom.c*u - geom.omega*time); ...
            geom.omega*geom.rf*cos(u).*sin(geom.c*u - geom.omega*time);];
    
        vel = v .* r0_t;
        
   end
   
   
