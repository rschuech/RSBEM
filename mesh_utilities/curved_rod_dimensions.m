function [minor_radius, major_radius, angle, arclength] = curved_rod_dimensions(AR1, AR2, V)

% computes dimensions of a sphere, straight rod (capsule) or curved rod
% (capped by hemispheres) based on AR1 = elongation (1 - Inf), AR2 = donutness (0 - < 1), V = volume

%AR1, AR2, V can be matrices of equal size, or if any of them is a scalar,
%it is assumed to be constant for all choices of the other parameters

%minor_radius = radius of sphere, cylinder, or torus section

%major_radius = radius of curvature, = Inf for spheres and straight rods,
%finite value for curved rods

%angle = NaN for spheres and straight rods, 0 - < 2*pi for curved rods (note
%that due to hemisphere caps, upper limit is actually < 2*pi before
%self-intersection occurs, but the code doesn't check for this)

%arclength = 0 for spheres, height of cylinder for straight rods, arc
%length of torus section centerline for curved rods (useful if making
%straight rods with AR2 = 0)


%% boring input expansions
if ~isscalar(AR1)
    siz = size(AR1);
elseif ~isscalar(AR2)
    siz = size(AR2);
else
    siz = size(V);
end

if isscalar(AR1)
    AR1 = repmat(AR1, siz);
end

if isscalar(AR2)
    AR2 = repmat(AR2, siz);
end

if isscalar(V)
    V = repmat(V, siz);
end
%%    


minor_radius = NaN(size(AR1));
major_radius = minor_radius;
angle = minor_radius;
arclength = minor_radius;



%OscarKatie E.coli is ~ 1.9 um long and 0.9 um wide, giving volume ~ 1
%um^3, so this is used for all bodies

sphere_radius = (V*3/4/pi).^(1/3);  %equivalent sphere radius


parfor cc = 1:numel(AR1)
    %radius = [];
 
    if AR1(cc) == 1  %we have a sphere, AR2 can only be 0
        
        if AR2(cc) == 0
            minor_radius(cc) = sphere_radius(cc);
            arclength(cc) = 0;
            angle(cc) = NaN;
            major_radius(cc) = Inf;
            
        else %you screwed up your inputs
            minor_radius(cc) = NaN;
            arclength(cc) = NaN;
            angle(cc) = NaN;
            major_radius(cc) = NaN;
        end
        
    elseif AR2(cc) == 0  % don't have a sphere, but do have a straight rod of zero curvature
        
        temp = roots([2/3*pi*(3*AR1(cc) - 1) 0 0 -V(cc)]); %see tablet worksheet for derivation
        for i = 1:length(temp)
            if isreal(temp(i)) && temp(i) > 0
                minor_radius(cc) = temp(i);  % major radius of cross section ellipse
                break
            end
        end
        
        if sum(imag(temp) == 0) > 1
            error('Problem with solving for geometry parameters');
        end
        
        arclength(cc) = 2*minor_radius(cc)*(AR1(cc) - 1);  %height of cylinder
        if arclength(cc) < 0
            error(['Calculations yield negative height = ',num2str(arclength(cc))]);
        end
      
        major_radius(cc) = Inf;  %this is how we inform Salome that this is a straight rod
        angle(cc) = NaN;  %meaningless for straight rod

    else % have a general curved rod
        
        temp = roots([2*pi*(AR1(cc) - 1/3) 0 0 -V(cc)]);
        for i = 1:length(temp)
            if isreal(temp(i)) && temp(i) > 0
                minor_radius(cc) = temp(i);
                break
            end
        end
        
        arclength(cc) = 2*minor_radius(cc)*(AR1(cc) - 1);  %arclength of torus section centerline
        major_radius(cc) = (arclength(cc) + 2*minor_radius(cc)) / AR2(cc) / 2 / pi;
        
        if major_radius(cc) - minor_radius(cc) < eps  %if radius of inner boundary is nearly zero or negative, this is a horn or spindle torus and self-intersecting
            arclength(cc) = NaN;
            minor_radius(cc) = NaN;
            major_radius(cc) = NaN;
            angle(cc) = NaN;
        end
        
        angle(cc) = arclength(cc) / (major_radius(cc));  %in radians
        
    end
    
end
    
    