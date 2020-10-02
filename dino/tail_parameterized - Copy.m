function [pt, vel, der, V_u, radius] = tail_parameterized(t, time, geom)
%just to be extra confusing, t is used as the parameter for position along
%tail and time is used for, duh, time.  geometric tail parameters stored in geom.

%tail is formed of two exponentially decaying sections continuously joined
%at t_transition

%optional 2nd output is velocity of each pt
%radius is local radius, including points on end spheres (starting and
%ending!!)

pt = NaN(length(t),3);  vel = pt;  der = pt;  V_u = pt;

radius = NaN(length(t),1);


for i = 1:length(t)
    
    
    %     dt = t(i) - geom.t_transition;
    %     d = diff(geom.amp);
    %
    %     C1 = 1/2*geom.kE(1)^2*geom.lambda^2*t(i)/pi^2;
    %     B1 = geom.amp(1)*C1*exp(-C1*t(i)/2);
    %     B2 = geom.amp(1)*(1-exp(-C1*t(i)/2));
    
    %         if t(i) > geom.t_transition
    %     C2 = 1/2*geom.kE(2)^2*geom.lambda^2*dt/pi^2;
    %     B3 = d*C2*exp(-C2*t(i)/2);
    %     B4 = d*(1-exp(-C2*t(i)/2));
    %         end
    if t(i) <= geom.t_min  %must be in starting sphere
        
        start_der = [ -geom.lambda/(2*pi)  0  0]';
        start_pt = [0 0 0]';  %appears to be the case according to eq for t(i) <= geom.t_transition with t = 0
        
        delta = ( t(i) - geom.t_min) *(start_der);  %what to add to start_pt to get to linearly extrapolated point inside starting sph
        pt(i,:) = start_pt + delta;  %parametric eq for a line "starting" at start_pt with slope start_der
        
        dist = min(geom.radius, sqrt(sum(delta.^2)));  %distance from start_pt that we must move to get to linearly extrapolated pt inside end sph
        radius(i) = sqrt( geom.radius^2 - dist^2);  %equation for local radius is basically eq of sphere with y = 0, giving z = +/- sqrt(r^2 - x^2) where x is distance from sph center i.e. end_pt
        
    elseif t(i) <= geom.t_transition
        pt(i,:) = [-geom.lambda*t(i)/(2*pi), ...
            -( geom.amp(1)*(1-exp(-(geom.kE(1))^2 * (geom.lambda*t(i)/2/pi)^2))*sin(t(i) - geom.omega*time) ), ...
            0];
        radius(i) = geom.radius;
    elseif t(i) > geom.t_transition && t(i) <= geom.t_max
        pt(i,:) =  [-geom.lambda*t(i)/(2*pi), ...
            -( ( geom.amp(1)*(1-exp(-(geom.kE(1))^2 * (geom.lambda*t(i)/2/pi)^2)) + diff(geom.amp)*(1 - exp(-(geom.kE(2))^2 * (geom.lambda*(t(i) - geom.t_transition)/2/pi)^2)) ) * sin(t(i) - geom.omega*time)  ) ,...
            0];
        radius(i) = geom.radius;
    else
        end_pt =  [-geom.lambda*geom.t_max/(2*pi), ...
            -( ( geom.amp(1)*(1-exp(-(geom.kE(1))^2 * (geom.lambda*geom.t_max/2/pi)^2)) + diff(geom.amp)*(1 - exp(-(geom.kE(2))^2 * (geom.lambda*(geom.t_max - geom.t_transition)/2/pi)^2)) ) * sin(geom.t_max - geom.omega*time)  ) ,...
            0];
        end_der =  [-geom.lambda/(2*pi), ...
            -(  -((1/2)*geom.amp(1)*geom.kE(1)^2*geom.lambda^2*geom.t_max*exp(-(1/4)*geom.kE(1)^2*geom.lambda^2*geom.t_max^2/pi^2)/pi^2+(1/2)*diff(geom.amp)*geom.kE(2)^2*geom.lambda^2*(geom.t_max - geom.t_transition)*...
            exp(-(1/4)*geom.kE(2)^2*geom.lambda^2*(geom.t_max- geom.t_transition)^2/pi^2)/pi^2)*sin(-geom.t_max+geom.omega*time)+...
            (geom.amp(1)*(1-exp(-(1/4)*geom.kE(1)^2*geom.lambda^2*geom.t_max^2/pi^2))+diff(geom.amp)*(1-exp(-(1/4)*geom.kE(2)^2*geom.lambda^2*(geom.t_max-geom.t_transition)^2/pi^2)))*cos(-geom.t_max+geom.omega*time)  ) ,...
            0];
        delta = ( t(i) - geom.t_max) *(end_der);  %what to add to end_pt to get to linearly extrapolated point inside end sph
        pt(i,:) = end_pt + delta;  %parametric eq for a line "starting" at end_pt with slope end_der
        % pt(i,:) = end_pt + t(i)*(end_der / sqrt(sum(end_der.^2)));  %parametric eq for a line "starting" at end_pt with slope end_der
        dist = min(geom.radius, sqrt(sum(delta.^2)));  %distance from end_pt that we must move to get to linearly extrapolated pt inside end sph
        radius(i) = sqrt( geom.radius^2 - dist^2);  %equation for local radius is basically eq of sphere with y = 0, giving z = +/- sqrt(r^2 - x^2) where x is distance from sph center i.e. end_pt
    end
    
    
    %rotation by 180 deg around z axis (to flip direction of tail along
    %x-axis) has already been incorporated into above eqs (negative x
    %component vs Salome)
    
    % at this point, tail is pointed in -x dir and is in x-y plane
    
    %now rotate around axis defined by point geom.rotation_pt and direction
    %geom.rotation_vec by angle geom.tail_angle  degrees
    
    %http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
    for ri = 1:size(geom.rotation_vec,2)   %currently, only consists of 90 deg around x axis to flip tail into x-z plane
        pt(i,:) = rotate_arbitrary_axis(pt(i,:), geom.rotation_pt(:,ri), geom.rotation_vec(:,ri), geom.rotation_angle(ri));
    end
    
    
    %     pt(i,:) = rotate_arbitrary_axis(pt(i,:), geom.rotation_pt(:,ri), geom.rotation_vec(:,ri), geom.tail_angle * pi/180);
    
    
    pt(i,:) = pt(i,:) + geom.translation';
    
    
    if nargout >= 2 %also get velocity
        
        if t(i) <= geom.t_min
            vel(i,:) = [0 0 0];
            
        elseif t(i) <= geom.t_transition
            vel(i,:) = [0 ,...
                -(  geom.amp(1)*(1-exp(-(geom.kE(1))^2 * (geom.lambda*t(i)/2/pi)^2))*(-geom.omega)*cos(t(i) - geom.omega*time) ) ,...
                0];
        elseif t(i) > geom.t_transition && t(i) <= geom.t_max
            vel(i,:) =  [0 ,...
                -(  ( geom.amp(1)*(1-exp(-(geom.kE(1))^2 * (geom.lambda*t(i)/2/pi)^2)) + diff(geom.amp)*(1 - exp(-(geom.kE(2))^2 * (geom.lambda*(t(i) - geom.t_transition)/2/pi)^2)) ) * (-geom.omega)*cos(t(i) - geom.omega*time)  ), ...
                0];
            
        else
            %               temp1 =  [0, ...
            %               geom.omega *  ( ( geom.amp(1)*(1-exp(-(geom.kE(1))^2 * (geom.lambda*geom.t_max/2/pi)^2)) + diff(geom.amp)*(1 - exp(-(geom.kE(2))^2 * (geom.lambda*(geom.t_max - geom.t_transition)/2/pi)^2)) ) * cos(geom.t_max - geom.omega*time)  ) ,...
            %             0];
            %
            %         tenp2 =  [0, ...
            %       ( t(i) - geom.t_max) *  ( geom.omega * (B1 + B3)*cos(-geom.t_max + geom.omega * time) + geom.omega * (B2 + B4)*sin(-geom.t_max + geom.omega * time)    );
            %                 0];
            dt = geom.t_max - geom.t_transition;
            d = diff(geom.amp);
            
            C1 = 1/2*geom.kE(1)^2*geom.lambda^2*geom.t_max/pi^2;
            B1 = geom.amp(1)*C1*exp(-C1*geom.t_max/2);
            B2 = geom.amp(1)*(1-exp(-C1*geom.t_max/2));
            
            
            C2 = 1/2*geom.kE(2)^2*geom.lambda^2*dt/pi^2;
            B3 = d*C2*exp(-C2*geom.t_max/2);
            B4 = d*(1-exp(-C2*geom.t_max/2));
            
            vel(i,:) = [0, ...
                geom.omega *  ( ( geom.amp(1)*(1-exp(-(geom.kE(1))^2 * (geom.lambda*geom.t_max/2/pi)^2)) + diff(geom.amp)*(1 - exp(-(geom.kE(2))^2 * (geom.lambda*(geom.t_max - geom.t_transition)/2/pi)^2)) ) * cos(geom.t_max - geom.omega*time)  ) + ...
                ( t(i) - geom.t_max) *  ( geom.omega * (B1 + B3)*cos(-geom.t_max + geom.omega * time) + geom.omega * (B2 + B4)*sin(-geom.t_max + geom.omega * time)    ), ...
                0];
        end
        
        % vel(i,:) = rotate_arbitrary_axis(vel(i,:), [0 0 0]', geom.rotation_vec, geom.tail_angle * pi/180);  %these are vectors, so make rotation point the origin
        
        for ri = 1:size(geom.rotation_vec,2)   %currently, only consists of 90 deg around x axis to flip tail into x-z plane
            vel(i,:) = rotate_arbitrary_axis(vel(i,:), geom.rotation_pt(:,ri), geom.rotation_vec(:,ri), geom.rotation_angle(ri));
        end
        
    end
    
    
    if nargout >= 3
        
        if t(i) <= geom.t_min
            der(i,:) = start_der;
        elseif t(i) <= geom.t_transition
            der(i,:) = [-geom.lambda/(2*pi), ...
                -( geom.amp(1)*( 1/2*geom.kE(1)^2*geom.lambda^2*t(i)/(pi^2)*exp(-geom.kE(1)^2*(geom.lambda*t(i)/2/pi)^2)*sin(t(i) - geom.omega*time) ...
                + (1-exp(-(geom.kE(1))^2 * (geom.lambda*t(i)/2/pi)^2))*cos(t(i) - geom.omega*time)) ), ...
                0];
            
            
        elseif t(i) > geom.t_transition && t(i) <= geom.t_max
            
            der(i,:) =  [-geom.lambda/(2*pi), ...
                -(  -((1/2)*geom.amp(1)*geom.kE(1)^2*geom.lambda^2*t(i)*exp(-(1/4)*geom.kE(1)^2*geom.lambda^2*t(i)^2/pi^2)/pi^2+(1/2)*diff(geom.amp)*geom.kE(2)^2*geom.lambda^2*(t(i) - geom.t_transition)*...
                exp(-(1/4)*geom.kE(2)^2*geom.lambda^2*(t(i)- geom.t_transition)^2/pi^2)/pi^2)*sin(-t(i)+geom.omega*time)+...
                (geom.amp(1)*(1-exp(-(1/4)*geom.kE(1)^2*geom.lambda^2*t(i)^2/pi^2))+diff(geom.amp)*(1-exp(-(1/4)*geom.kE(2)^2*geom.lambda^2*(t(i)-geom.t_transition)^2/pi^2)))*cos(-t(i)+geom.omega*time)  ) ,...
                0];
        else
            
            %   pt(i,:) = end_pt + ( t(i) - geom.t_max) *(end_der);
            %   derivative WRT t(i) is simply end_der (a constant, no more
            %   dependence on t, but it still depends on time)
            
            der(i,:) = end_der;  %should have been computed above
            
        end
        
        %der(i,:) = rotate_arbitrary_axis(der(i,:), [0 0 0]', geom.rotation_vec, geom.tail_angle * pi/180);  %these are vectors, so make rotation point the origin
        
        for ri = 1:size(geom.rotation_vec,2)   %currently, only consists of 90 deg around x axis to flip tail into x-z plane
            der(i,:) = rotate_arbitrary_axis(der(i,:), geom.rotation_pt(:,ri), geom.rotation_vec(:,ri), geom.rotation_angle(ri));
        end
        
        % "upward" pointing vector = rotation of tangent der vector by -90 deg
        % around z-axis
        %
       % V_u(i,:) = [der(i,2), -der(i,1), der(i,3)];
                % "upward" pointing vector = rotation of tangent der vector by 90 deg
        % around y-axis
        %
            V_u(i,:) = [der(i,3), -der(i,2), -der(i,1)];
    end
    
    
    
    %     if t(i) <= geom.t_transition
    %         V_u(i,:) = [- ( geom.amp(1)*( 1/2*geom.kE(1)^2*geom.lambda^2*t(i)/(pi^2)*exp(-geom.kE(1)^2*(geom.lambda*t(i)/2/pi)^2)*sin(t(i) - geom.omega*time) ...
    %             + (1-exp(-(geom.kE(1))^2 * (geom.lambda*t(i)/2/pi)^2))*cos(t(i) - geom.omega*time)) ) ,...
    %             geom.lambda/(2*pi),...
    %             0];
    %
    %
    %
    %   elseif t(i) > geom.t_transition && t(i) <= geom.t_max
    %
    %         V_u(i,:) =  [ -(  -((1/2)*geom.amp(1)*geom.kE(1)^2*geom.lambda^2*t(i)*exp(-(1/4)*geom.kE(1)^2*geom.lambda^2*t(i)^2/pi^2)/pi^2+(1/2)*diff(geom.amp)*geom.kE(2)^2*geom.lambda^2*(t(i) - geom.t_transition)*...
    %             exp(-(1/4)*geom.kE(2)^2*geom.lambda^2*(t(i)- geom.t_transition)^2/pi^2)/pi^2)*sin(-t(i)+geom.omega*time)+...
    %             (geom.amp(1)*(1-exp(-(1/4)*geom.kE(1)^2*geom.lambda^2*t(i)^2/pi^2))+diff(geom.amp)*(1-exp(-(1/4)*geom.kE(2)^2*geom.lambda^2*(t(i)-geom.t_transition)^2/pi^2)))*cos(-t(i)+geom.omega*time)  ) ,...
    %             geom.lambda/(2*pi), ...
    %             0];
    %
    %
    %
    %     else
    %     V_u(i,:) = [end_der(2), -end_der(1), end_der(3)];
    %     end
    
    %V_u(i,:) = [der(i,2), -der(i,1), der(i,3)];
    
    
    
end

