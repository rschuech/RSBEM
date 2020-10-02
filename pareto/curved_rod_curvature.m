function [curvs] = curved_rod_curvature(AR1, AR2, V)
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



avg_mean_curv = NaN(size(AR1));
avg_abs_mean_curv = avg_mean_curv;
avg_abs_mean_abs_curv = avg_mean_curv;

avg_gauss_curv = avg_mean_curv;
avg_abs_gauss_curv = avg_mean_curv;

avg_max_abs_principle_curv = avg_mean_curv;
max_principle_curv = avg_mean_curv;
max_mean_curv = avg_mean_curv;
max_gauss_curv = avg_mean_curv;



%OscarKatie E.coli is ~ 1.9 um long and 0.9 um wide, giving volume ~ 1 um^3
% V = 1; %all body shapes will have this volume = 1 um^3
sphererad = (V*3/4/pi)^(1/3);  %equivalent sphere radius

ppm = ParforProgressStarter2('computing curvature metrics', numel(AR1));
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
    
    % below max values work for both curved and straight since radius_curv
    % = Inf describes straight
    
    % after much thought I can only guess that k1 an k2 of a torus must be
    % +1/a and +cos(v)/(c+a*cos(v) because the 1/a must be positive
    % curvature.  The sign on the other one doesn't really matter since it
    % depends on where the v coord system starts which is never defined,
    % but seems most straightforward to have it postive as well
    % this site seems to confirm:  http://physics.oregonstate.edu/coursewikis/GDF/book/gdf/torus
    
    max_principle_curv(cc) = max(1/radius, abs(1/(radius_curv - radius)));  %surface maximum of abs principle curvatures
    
    max_gauss_curv(cc) = max( 1/(radius*(radius_curv - radius))  ,   1/radius^2);   %surface max of abs gaussian curvature
    % max principle curv k1 for torus = 1/(radius_curv - radius) so max
    % gaussian curv for torus = k1*k2 = 1/(radius_curv - radius) * 1/radius
    
    abs_max_mean_curv_torus = max( abs([ (2*radius - radius_curv)/2/radius/(radius-radius_curv)      (2*radius+radius_curv)/2/radius/(radius+radius_curv)  ]) );  % see mean_curvature_torus.m - the extreme values of mean curv always occur at one of these two limits
    max_mean_curv(cc) = max(  abs_max_mean_curv_torus  ,  1/radius  );  %surface max of mean curvature -
    % abs max of mean curv over torus is the abs max mean curv of torus,
    % then over entire surface take the max of the
    % torus mean curv and the mean curv of the end spheres, 1/a
    
%     A = (4*pi*radius^2  +  (l-2*radius)*2*pi*radius); % total area - use
%     this to do surface averages of curvature
    A = 1;  % just do surface integrals of curvature, not averages - should really rename all the below variables since they're not avg in this case....
    
    if ~isinf(radius_curv) && ~isnan(radius_curv)  %curved rod
        %         mean_curv(cc) = 2/l + 1/2/l*( (2*sqrt(radius_curv^2 - radius^2) - radius_curv ) / (radius*sqrt(radius_curv^2 - radius^2)) ) *(l-2*radius); %avg value of mean curv over surface, allowing regions of negative curvature to offset regions of positive curvature
        
        fun = @(v) (   (radius_curv + 2.*radius.*cos(v))./(2.*radius.*(radius_curv+radius.*cos(v)))    );
        I = ( integral(fun,0,2*pi,'abstol',1E-12,'reltol',1E-10) );  %just in case integral is negative, take abs value
        %         avg_mean_curv(cc) = abs(  2/l + I*(l-2*radius)/2/pi/l   );  %surface average of regular old mean curv, abs taken only at very end to give postive cost value (I might be negative and cancel with 1/a for spheres)
        avg_mean_curv(cc) = abs( 1/radius * 4*pi*radius^2     +     I/(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
        fun = @(v) abs(   (radius_curv+2.*radius.*cos(v))./(2.*radius.*(radius_curv+radius.*cos(v)))     );  %avg value of abs(mean curv) so both neg and pos mean curvature add to the surface average, but at a single point, oppositely signed principle curvatures can still cancel
        I = abs(  integral(fun,0,2*pi,'abstol',1E-12,'reltol',1E-10)  );
        %         avg_abs_mean_curv(cc) = abs(  2/l + I*(l-2*radius)/2/pi/l  );  %surface average of absolute value of mean curv
        avg_abs_mean_curv(cc) = abs( 1/radius * 4*pi*radius^2     +     I/(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
        fun = @(v) abs(  1/2*( abs(1./radius) + abs(cos(v)./(radius_curv + radius.*cos(v))) ) );  %each principle curvature is absolute valued, so no negative values can occur anywhere
        I = abs(  integral(fun,0,2*pi,'abstol',1E-12,'reltol',1E-10) );
        %         avg_abs_mean_abs_curv(cc) = abs(  2/l + I*(l-2*radius)/2/pi/l  ); %surface average of absolute value (redundant) of mean abs curv (modified mean curv)
        avg_abs_mean_abs_curv(cc) = abs( 1/radius * 4*pi*radius^2     +     I/(2*pi) * (l-2*radius)*2*pi*radius    ) /A;
        
        fun = @(v) max( abs(1./radius) , abs(cos(v)./(radius_curv + radius.*cos(v))) );  %each principle curvature is absolute valued, and max taken at each point, then surface averaged
        I = abs( integral(fun,0,2*pi,'abstol',1E-12,'reltol',1E-10) );
        %         avg_max_abs_principle_curv(cc) = abs(  2/l + I*(l-2*radius)/2/pi/l  );
        avg_max_abs_principle_curv(cc) = abs( 1/radius * 4*pi*radius^2     +     I/(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
        fun = @(v) abs(    cos(v)./(radius_curv + radius.*cos(v))   .*    1/radius    );
        I = abs( integral(fun,0,2*pi,'abstol',1E-12,'reltol',1E-10) );
        %         avg_abs_gauss_curv(cc) = (4*pi + I*(l-2*radius)*radius) / (4*pi*radius^2 + (l-2*radius)*2*pi*radius);  %surface avg of abs value of gaussian curvature
        %          avg_abs_gauss_curv(cc) = (1/radius^2 * 4*pi*radius^2     +       I *(l-2*radius)*radius) ;  % gauss curv integrated over surface; unitless
        avg_abs_gauss_curv(cc) = abs( 1/radius^2 * 4*pi*radius^2     +     I/(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
        fun = @(v) (    cos(v)./(radius_curv + radius.*cos(v))   .*    1/radius    );
        I = (  integral(fun,0,2*pi,'abstol',1E-12,'reltol',1E-10)  );
        %         avg_gauss_curv(cc) = abs(   (4*pi + I*(l-2*radius)*radius) / (4*pi*radius^2 + (l-2*radius)*2*pi*radius)  );  %surface avg of gaussian curvature - abs only taken at very end to ensure positive cost
        avg_gauss_curv(cc) = abs( 1/radius^2 * 4*pi*radius^2     +     I/(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
        
    else  %straight rod
        %         mean_curv(cc) = 2/l + 1/2/l*(1/radius)*(l-2*radius); %shatlab doesn't know how to take the limit
        
        
        %         avg_mean_curv(cc) = 2/l + (l-2*radius)/2/radius/l;
        avg_mean_curv(cc) = abs( 1/radius * 4*pi*radius^2     +     1/(2*radius) *2*pi /(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
        %         avg_abs_mean_curv(cc) = 2/l + (l-2*radius)/2/radius/l;   %same as avg_mean_curv because in the limit of radius_curv = Inf, cos(v) doesn't matter
        avg_abs_mean_curv(cc) = abs( 1/radius * 4*pi*radius^2     +     1/(2*radius)*2*pi  /(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
        %         avg_abs_mean_abs_curv(cc) = 2/l + (l-2*radius)/2/radius/l; % again the same because in the limit of radius_curv = Inf, doesn't matter whether abs is taken or not, shat is always positive
        avg_abs_mean_abs_curv(cc) = abs( 1/radius * 4*pi*radius^2     +     1/(2*radius)*2*pi  /(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
        
        
        %         avg_max_abs_principle_curv(cc) = abs(  2/l + 1/radius *(l-2*radius)/l  );  % max abs princ curv on cylinder is 1/a, so I = 2*pi/a in eq for curved rod above
        avg_max_abs_principle_curv(cc) =  abs( 1/radius * 4*pi*radius^2     +     1/(radius)*2*pi  /(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
        %          avg_abs_gauss_curv(cc) = 2  / (2*radius^2 + (l-2*radius)*radius);  % gauss curv on cylinder = 0 so I = 0 in curved rod eq above
        avg_abs_gauss_curv(cc) = abs( 1/radius^2 * 4*pi*radius^2     +     0/(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
        
        %             avg_gauss_curv(cc) = 2  / (2*radius^2 + (l-2*radius)*radius);  % same as avg_abs_gauss_curv, since I = 0 still
        avg_gauss_curv(cc) = abs( 1/radius^2 * 4*pi*radius^2     +     0/(2*pi) * (l-2*radius)*2*pi*radius    ) / A;
        
    end
    ppm.increment(cc);
end
delete(ppm);

curvs.avg_mean_curv = avg_mean_curv;
curvs.avg_abs_mean_curv = avg_abs_mean_curv;
curvs.avg_abs_mean_abs_curv = avg_abs_mean_abs_curv;

curvs.avg_gauss_curv = avg_gauss_curv;
curvs.avg_abs_gauss_curv = avg_abs_gauss_curv;

curvs.avg_max_abs_principle_curv = avg_max_abs_principle_curv;
curvs.max_principle_curv = max_principle_curv;
curvs.max_mean_curv = max_mean_curv;
curvs.max_gauss_curv = max_gauss_curv;
