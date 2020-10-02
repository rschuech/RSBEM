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
   
%    r0 = [  ( geom.b*(u) + geom.rf * sin(geom.c*u - geom.omega*time) ) ; ...
%           (geom.R+geom.d+geom.w+geom.rf)*sin(u) + geom.rf* sin(u).*cos(geom.c*u - geom.omega*time); ...
%           -( (geom.R+geom.d+geom.w+geom.rf)*cos(u) +  geom.rf*cos(u).*cos(geom.c*u - geom.omega*time) );];
% 
%    pt = v.*r0 + [ (1-v).*geom.b.*u;  (geom.R+geom.d).*(1-v).*sin(u);   -( (geom.R+geom.d)*(1-v).*cos(u) );];
%    
%    pt(1) = pt(1) + geom.shift; %translate in x dir
%    
%    if nargout > 1
%         r0_t = [  ( - geom.omega*geom.rf * cos(geom.c*u - geom.omega*time) ) ; ...
%             geom.omega*geom.rf* sin(u).*sin(geom.c*u - geom.omega*time); ...
%             -( geom.omega*geom.rf*cos(u).*sin(geom.c*u - geom.omega*time) );];
%     
%         vel = v .* r0_t;
%         
%    end
   
   
   
   
% to be input into Salome, with spiral
% x = v * ( (R+d+w+rf)*cos(t) + rf* cos(t)*cos(c*t - omega*time) ) + (R+d)*(1-v).*cos(t) % (v, transverse_R, transverse_d, transverse_w, transverse_rf, transverse_rf, transverse_c, transverse_omega, time, transverse_R, transverse_d, v)
% y = v * ( (R+d+w+rf)*sin(t) + rf* sin(t)*cos(c*t - omega*time) ) + (R+d)*(1-v).*sin(t) % (v, transverse_R, transverse_d, transverse_w, transverse_rf, transverse_rf, transverse_c, transverse_omega, time, transverse_R, transverse_d, v)
% z = v * ( -b*(t) +             rf*        sin(c*t - omega*time) ) +     -(1-v).*b.*t % (v, transverse_b, transverse_rf, transverse_c, transverse_omega, time, v, transverse_b)

   
   %theta = 45 * pi/180;  %now this is again rotated by some random angle around x-axis as another Salome workaround....
   theta = 0;  %apparently don't need to rotate transverse for normal groove, only for big groove...
   
   
      r0 = [  ( geom.b*(u) + geom.rf * sin(geom.c*u - geom.omega*time) ) ; ...
          (geom.R+geom.d+geom.w+geom.rf)*sin(u) + geom.rf* sin(u).*cos(geom.c*u - geom.omega*time); ...
          -( (geom.R+geom.d+geom.w+geom.rf)*cos(u) +  geom.rf*cos(u).*cos(geom.c*u - geom.omega*time) );];

      
      rotmat = [1 0 0; 0 cos(theta) -sin(theta);  0 sin(theta) cos(theta)];
      
      r0 = rotmat * r0;
      temp = rotmat * [ (1-v).*geom.b.*u;  (geom.R+geom.d).*(1-v).*sin(u);   -( (geom.R+geom.d)*(1-v).*cos(u) );];
      
   pt = v.*r0 + temp;
   
   pt(1) = pt(1) + geom.shift; %translate in x dir
   
   if nargout > 1
        r0_t = [  ( - geom.omega*geom.rf * cos(geom.c*u - geom.omega*time) ) ; ...
            geom.omega*geom.rf* sin(u).*sin(geom.c*u - geom.omega*time); ...
            -( geom.omega*geom.rf*cos(u).*sin(geom.c*u - geom.omega*time) );];
    
        r0_t = rotmat * r0_t;
        
        vel = v .* r0_t;
        
   end
