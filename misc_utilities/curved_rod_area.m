function [area, radius_curv1, radius_curv2] = curved_rod_area(AR1, AR2, V)


area = NaN(size(AR1));
radius_curv1 = area;
radius_curv2 = area;



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
        
        area(cc) = 4*pi*radius^2;
        
        radius_curv1(cc) = radius;
        radius_curv2(cc) = radius_curv;
        
        % cost = curvature * area
        % ease = 1/curv * 1/area     =     radius_curv / area
%         cyl_area = 0;
%         radius_curv1(cc) = radius / cyl_area;
%         radius_curv2(cc) = Inf;
        
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
        
        area(cc) = 4*pi*radius^2 + 2*pi*radius*height;
        radius_curv1(cc) = radius;
        radius_curv2(cc) = Inf;
%         
%         cyl_area = 2*pi*radius*height;
%         radius_curv1(cc) = radius / cyl_area;
%         radius_curv2(cc) = Inf;
        
    else % have a general curved rod
        
        temp = roots([2*pi*(AR1(cc) - 1/3) 0 0 -V]);
        for i = 1:length(temp)
            if isreal(temp(i)) && temp(i) > 0
                radius = temp(i);
                break
            end
        end
        
        arclength = 2*radius*AR1(cc);  % total arclength to ends of spheres
        radius_curv = arclength / AR2(cc) / 2 / pi;
        
        if radius_curv - radius < eps  %if radius of inner boundary is zero or negative, this shape is impossible
            arclength = NaN;
            radius_curv = NaN;
            radius = NaN;
        end
        
        
        nturns = (arclength - 2 * radius) / 2 / pi / radius_curv;  %will always be a fraction - using modified helix-creating code to generate curved rods in a plane
        height = NaN;
        
        l = arclength;
        
        area(cc) = 4*pi*radius^2 + 2*pi*radius*(arclength - 2*radius);
        radius_curv1(cc) = radius;
        radius_curv2(cc) = radius_curv;
        
%         cyl_area = 2*pi*radius*(arclength - 2*radius);
%         radius_curv1(cc) = radius / cyl_area;
%         radius_curv2(cc) = radius_curv / cyl_area;
%         
    end
    
    
end




