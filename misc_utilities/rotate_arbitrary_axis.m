function rotated = rotate_arbitrary_axis(input, point, vec, angle)

    vec = vec / sqrt(sum(vec.^2));  %make sure it's a unit vec
    
    a = point(1);  b = point(2);  c = point(3);
    u = vec(1);  v = vec(2);  w = vec(3);
    
    
temp = input; 
parfor i = 1:size(temp,1)
    input = temp(i,:);

    
    rotated(i,:) = [ (a*(v^2+w^2) - u*(b*v+c*w-u*input(1)-v*input(2)-w*input(3)))*(1-cos(angle)) + input(1)*cos(angle) + (-c*v + b*w - w*input(2) + v*input(3)) * sin(angle); ...
        (b*(u^2+w^2) - v*(a*u+c*w-u*input(1)-v*input(2)-w*input(3)))*(1-cos(angle)) + input(2)*cos(angle) + (c*u - a*w + w*input(1) - u*input(3)) * sin(angle); ...
        (c*(u^2+v^2) - w*(a*u+b*v-u*input(1)-v*input(2)-w*input(3)))*(1-cos(angle)) + input(3)*cos(angle) + (-b*u + a*v - v*input(1) + u*input(2))*sin(angle); ];
    
end