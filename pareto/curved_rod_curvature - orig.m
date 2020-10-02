function [mean_curv, max_curv, mean_abs_curv, mean_max_curv] = curved_rod_curvature(AR1, AR2, V)
% computes mean and max curvature for curved rods, straight rods, and
% spheres

% mean curvature is the average value of the mean curvature over the
% surface

% max curvature is the absolute max principle curvature anywhere on the
% surface

% mean abs curvature is the average value over the surface of the mean of the absolute
% principle curvatures.  this is different from mean curvature because for
% a torus, the toroidal principle curvature component is mostly negative, so can
% cancel out the cylindrical principle curvature component.  Taking
% absolute values before the mean results in a modified mean curvature that
% doesn't allow any offsetting.

% mean max curvature is the average value over the surface of the absolute
% max principle curvature.



mean_curv = NaN(size(AR1));
max_curv = mean_curv;
mean_abs_curv = mean_curv; 
mean_max_curv = mean_curv;



%OscarKatie E.coli is ~ 1.9 um long and 0.9 um wide, giving volume ~ 1 um^3
% V = 1; %all body shapes will have this volume = 1 um^3
sphererad = (V*3/4/pi)^(1/3);  %equivalent sphere radius


parfor cc = 1:numel(AR1)
    radius = [];
     %cc/numel(AR1)
    if AR1(cc) == 1  %we have a sphere
        
        if AR2(cc) == 0
            radius = sphererad;
            height = 0;
            totlength = 2*sphererad;
            nturns = NaN;
            radius_curv = Inf;
            
            l = totlength;
            
        else
            radius = NaN;
            height = NaN;
            totlength = NaN;
            radius_curv = NaN;
            l = NaN;
        end
        
    elseif AR2(cc) == 0  % don't have a sphere, but do have a straight rod of zero curvature
        
        temp = roots([2/3*pi*(3*AR1(cc) - 1) 0 0 -V]); %see tablet worksheet
        for i = 1:length(temp)
            if isreal(temp(i)) && temp(i) > 0
                radius = temp(i);  % major radius of cross section ellipse
                break
            end
        end
        
        if sum(imag(temp) == 0) > 1
            error('Problem with solving for geometry parameters');
        end
        
        height = 2*radius*(AR1(cc) - 1);  %height of cylinder
        if height < 0
            error(['Calculations yield negative height = ',num2str(height)]);
        end
        totlength = height + 2*radius;  %length of complete capsule
        
        radius_curv = Inf;  %this is how we inform Salome that this is a straight rod
        nturns = NaN;
        
        l = totlength;
    else % have a general curved rod
        
        temp = roots([2*pi*(AR1(cc) - 1/3) 0 0 -V]);
        for i = 1:length(temp)
            if isreal(temp(i)) && temp(i) > 0
                radius = temp(i);
                break
            end
        end
        
        arclength = 2*radius*AR1(cc);
        radius_curv = arclength / AR2(cc) / 2 / pi;
        
        if radius_curv - radius < eps  %if radius of inner boundary is zero or negative, this shape is impossible
            arclength = NaN;
            radius_curv = NaN;
            radius = NaN;
        end
        
        
        nturns = (arclength - 2 * radius) / 2 / pi / radius_curv;  %will always be a fraction - using modified helix-creating code to generate curved rods in a plane
        height = NaN;
        
        l = arclength;
    end
    
    
    
    max_curv(cc) = max(1/radius, 1/(radius_curv - radius));
    
    if ~isinf(radius_curv) && ~isnan(radius_curv)
        %         mean_curv(cc) = 2/l + 1/2/l*( (2*sqrt(radius_curv^2 - radius^2) - radius_curv ) / (radius*sqrt(radius_curv^2 - radius^2)) ) *(l-2*radius); %avg value of mean curv over surface, allowing regions of negative curvature to offset regions of positive curvature
        
        
        fun = @(v) abs((radius_curv+2*radius*cos(v))./(2*radius*(radius_curv+radius*cos(v))));  %avg value of abs(mean curv) so both neg and pos mean curvature contribute to the surface average, but at a single point, oppositely signed principle curvatures can still cancel
        I = integral(fun,0,2*pi,'abstol',1E-12,'reltol',1E-10);
        
%         I = integral(fun,0,2*pi);
        
        mean_curv(cc) = 2/l + I*(l-2*radius)/2/pi/l;
        
        fun = @(v) 1/2*( abs(1./radius) + abs(cos(v)./(radius_curv + radius.*cos(v))) );  %each principle curvature is absolute valued, so no negative values can occur anywhere
        I = integral(fun,0,2*pi,'abstol',1E-12,'reltol',1E-10);
        
        mean_abs_curv(cc) = 2/l + I*(l-2*radius)/2/pi/l;
        
          fun = @(v) max( abs(1./radius) , abs(cos(v)./(radius_curv + radius.*cos(v))) );  %each principle curvature is absolute valued, and max taken at each point, then surface averaged
        I = integral(fun,0,2*pi,'abstol',1E-12,'reltol',1E-10);
        
        mean_max_curv(cc) = 2/l + I*(l-2*radius)/2/pi/l;
        
    else
        %         mean_curv(cc) = 2/l + 1/2/l*(1/radius)*(l-2*radius); %shatlab doesn't know how to take the limit
        
        mean_curv(cc) = 2/l + (l-2*radius)/2/radius/l;
        
           mean_abs_curv(cc) = 2/l + 1/radius*(l-2*radius)/2/l;
           
         mean_max_curv(cc) = 2/l + 1/radius*(l-2*radius)/l;
           
        
    end
    
end
