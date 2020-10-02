function [pt, vel, der, V_u] = tail_parameterized(t, time, geom)
%just to be extra confusing, t is used as the parameter for position along
%tail and time is used for, duh time.  geometric tail parameters stored in geom.

%tail is formed of two exponentially decaying sections continuously joined
%at t_transition

%optional 2nd output is velocity of each pt

pt = NaN(length(t),3);  vel = pt;  der = pt;  V_u = pt;



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
    
    if t(i) <= geom.t_transition
        pt(i,:) = [-geom.lambda*t(i)/(2*pi), ...
            -( geom.amp(1)*(1-exp(-(geom.kE(1))^2 * (geom.lambda*t(i)/2/pi)^2))*sin(t(i) - geom.omega*time) ), ...
            0];
    elseif t(i) > geom.t_transition && t(i) <= geom.t_max
        pt(i,:) =  [-geom.lambda*t(i)/(2*pi), ...
            -( ( geom.amp(1)*(1-exp(-(geom.kE(1))^2 * (geom.lambda*t(i)/2/pi)^2)) + diff(geom.amp)*(1 - exp(-(geom.kE(2))^2 * (geom.lambda*(t(i) - geom.t_transition)/2/pi)^2)) ) * sin(t(i) - geom.omega*time)  ) ,...
            0];
    else
        end_pt =  [-geom.lambda*geom.t_max/(2*pi), ...
            -( ( geom.amp(1)*(1-exp(-(geom.kE(1))^2 * (geom.lambda*geom.t_max/2/pi)^2)) + diff(geom.amp)*(1 - exp(-(geom.kE(2))^2 * (geom.lambda*(geom.t_max - geom.t_transition)/2/pi)^2)) ) * sin(geom.t_max - geom.omega*time)  ) ,...
            0];
        end_der =  [-geom.lambda/(2*pi), ...
            -(  -((1/2)*geom.amp(1)*geom.kE(1)^2*geom.lambda^2*geom.t_max*exp(-(1/4)*geom.kE(1)^2*geom.lambda^2*geom.t_max^2/pi^2)/pi^2+(1/2)*diff(geom.amp)*geom.kE(2)^2*geom.lambda^2*(geom.t_max - geom.t_transition)*...
            exp(-(1/4)*geom.kE(2)^2*geom.lambda^2*(geom.t_max- geom.t_transition)^2/pi^2)/pi^2)*sin(-geom.t_max+geom.omega*time)+...
            (geom.amp(1)*(1-exp(-(1/4)*geom.kE(1)^2*geom.lambda^2*geom.t_max^2/pi^2))+diff(geom.amp)*(1-exp(-(1/4)*geom.kE(2)^2*geom.lambda^2*(geom.t_max-geom.t_transition)^2/pi^2)))*cos(-geom.t_max+geom.omega*time)  ) ,...
            0];
        pt(i,:) = end_pt + ( t(i) - geom.t_max) *(end_der);  %parametric eq for a line "starting" at end_pt with slope end_der
        % pt(i,:) = end_pt + t(i)*(end_der / sqrt(sum(end_der.^2)));  %parametric eq for a line "starting" at end_pt with slope end_der
    end
    pt(i,:) = pt(i,:) + geom.translation';
    
    
    
    %now rotate around axis defined by point geom.rotation_pt and direction
    %geom.rotation_vec by angle geom.tail_angle  degrees
    
    %http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
    pt(i,:) = rotate_arbitrary_axis(pt(i,:), geom.rotation_pt, geom.rotation_vec, geom.tail_angle * pi/180);
    
    
    
    if nargout >= 2 %also get velocity
        if t(i) <= geom.t_transition
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
        
        vel(i,:) = rotate_arbitrary_axis(vel(i,:), [0 0 0]', geom.rotation_vec, geom.tail_angle * pi/180);  %these are vectors, so make rotation point the origin
        
    end
    
    
    if nargout >= 3
        
        
        if t(i) <= geom.t_transition
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
        
         der(i,:) = rotate_arbitrary_axis(der(i,:), [0 0 0]', geom.rotation_vec, geom.tail_angle * pi/180);  %these are vectors, so make rotation point the origin
        
        
            % "upward" pointing vector = rotation of tangent der vector by -90 deg
    % around z-axis
    %
      V_u(i,:) = [der(i,2), -der(i,1), der(i,3)];
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

