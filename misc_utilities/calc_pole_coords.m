function [poles, sphere_centers, centerline_center] = calc_pole_coords(Metadata)

% poles are the points on the body surface at the extreme ends of the
% centerline
% sphere_centers are the centers of the head and tail capping spheres
% centerline_center is the center of the circle traced by the body
% centerline

switch Metadata(1).geom.shape
    case 'curved_rod'
        if Metadata(1).geom.AR1 == 1 %sphere degenerate case, origin is body center
            % pole_separation = Metadata(1).geom.radius * 2;
            tail_pole_x = - Metadata(1).geom.radius;  tail_pole_y = 0;
            head_pole_x =   Metadata(1).geom.radius;  head_pole_y = 0;
            sphere_centers.head = [0 0 0]';  sphere_centers.tail = [0 0 0]';
            centerline_center = [0 Inf 0]';
        elseif Metadata(1).geom.AR2 == 0  %not a sphere, but a straight rod degenerate case, origin is body center
            % pole_separation = Metadata(1).geom.radius * 2 + Metadata(1).geom.height;
            tail_pole_x = - (Metadata(1).geom.height/2 + Metadata(1).geom.radius);  tail_pole_y = 0;
            head_pole_x =   (Metadata(1).geom.height/2 + Metadata(1).geom.radius);    head_pole_y = 0;
             sphere_centers.head = [Metadata(1).geom.height/2 0 0]';  sphere_centers.tail = [-Metadata(1).geom.height/2 0 0]';
             centerline_center = [0 Inf 0]';
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
            
              sphere_centers.head = [x_1 + a, y_1 + b, 0]';  sphere_centers.tail = [0 0 0]';
              centerline_center = [x_1 y_1 0]';
        end
        
    case 'ellipsoid'
        head_pole_x = Metadata(1).geom.a;
        head_pole_y = 0;
        tail_pole_x = - Metadata(1).geom.a;
        tail_pole_y = 0;
end

poles.head = [head_pole_x; head_pole_y; 0];
poles.tail = [tail_pole_x; tail_pole_y; 0];