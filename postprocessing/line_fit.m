function [rmse, fitted, dists] = line_fit(intercept, slope, speed, t, points)

x = intercept(1) + slope(1) * t * speed;
y = intercept(2) + slope(2) * t * speed;
z = intercept(3) + slope(3) * t * speed;

    fitted = [x y z];
    
    
    dists = sqrt( sum( ( fitted - points ).^2 , 2));
    rmse = sqrt(sum(dists.^2))  /  length(t) ;
