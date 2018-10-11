function [pole_separation] = calc_pole_separation(input, Metadata, avg_swimming_axis)

switch input.body.shape
    case 'curved_rod'
        if input.body.AR(1) == 1 %sphere degenerate case, origin is body center
            % pole_separation = Metadata(1).geom.radius * 2;
            tail_pole_x = - Metadata(1).geom.radius;  tail_pole_y = 0;
            head_pole_x =   Metadata(1).geom.radius;    head_pole_y = 0;
        elseif input.body.AR(2) == 0  %not a sphere, but a straight rod degenerate case, origin is body center
            % pole_separation = Metadata(1).geom.radius * 2 + Metadata(1).geom.height;
            tail_pole_x = - (Metadata(1).geom.height/2 + Metadata(1).geom.radius);  tail_pole_y = 0;
            head_pole_x =   (Metadata(1).geom.height/2 + Metadata(1).geom.radius);    head_pole_y = 0;
        else  %actually a general curved rod, origin is center of tail-end sphere
            
            theta = Metadata(1).geom.nturns * 2 * pi;
            a = Metadata(1).geom.radius_curv * cos( theta - pi/2);
            b = Metadata(1).geom.radius_curv * sin( theta - pi/2);
            r = Metadata(1).geom.radius;
            x_1 = 0; y_1 = Metadata(1).geom.radius_curv;
            
            head_pole_x = x_1 + a + r*cos(theta);
            head_pole_y = y_1 + b + r*sin(theta);
            
            tail_pole_x = - Metadata(1).geom.radius;
            tail_pole_y = 0;
            
        end
        
    case 'ellipsoid'
        head_pole_x = Metadata(1).geom.a;
        head_pole_y = 0;
        tail_pole_x = - Metadata(1).geom.a;
        tail_pole_y = 0;
end


%pole_separation = sqrt( (head_pole_x - tail_pole_x)^2 + (head_pole_y - tail_pole_y)^2 );
% actually want distance between projection of
% poles onto swimming path for fore-aft separation

slope = avg_swimming_axis;  % only thing that matters, we only need a direction vector, not the full parametric line

[pole_separation] = line_projection_distance([tail_pole_x tail_pole_y 0]', [head_pole_x head_pole_y 0]',slope);