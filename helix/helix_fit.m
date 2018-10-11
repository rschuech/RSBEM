function [rmse, fitted] = helix_fit(a1, a2, v, f, phase,shift, angles,t, points)

% a1 = amp 1 (elliptical helix) 
% a2 = amp 2 (elliptical helix)
% v  = speed along x-axis (on avg)
% f = frequency
% phase = phase shift
% shift = x y z translation
% angles = x y z rotation
% t = time data of data to fit
% points = x y z location data to fit

% helix centerline assumed to more or less coincide with x-axis


    
    % t = linspace(0,30,100);
    %a1 = 1;  a2 = 1;  v = 1;
    %f1 = 1;  f2 = 1;
    % phase = 0;
    
    x = v*t;
    y = a1*sin(f*t+phase);
    z = a2*cos(f*t+phase);
    
    rotmat = rotation_matrix('z',angles(3)) * rotation_matrix('y',angles(2)) * rotation_matrix('x',angles(1));
    rotated = [x y z] * rotmat;
    x = rotated(:,1);  y = rotated(:,2);  z = rotated(:,3);
    
    %shift = [0 1 0]';
    x = x + shift(1);
    y = y + shift(2);
    z = z + shift(3);
    
    
    %angles = [0 0 0];
    
    fitted = [x y z];
    
    
    dists = sqrt( sum( ( fitted - points ).^2 , 2));
    rmse = sqrt(sum(dists.^2));


% figure(34)
% plot3(x,y,z,'o-');  axis equal