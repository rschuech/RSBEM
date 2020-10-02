function [obj] = objfun_contours(coord1, coord2, tangent1, tangent2)

if any(isnan( [ coord1(:); coord2(:); tangent1; tangent2]))
    obj = NaN;
else
    
    dist = sqrt( sum( ( coord1 - coord2).^2) );
    
    
%     derdiff = abs(der1 - der2);

dotprod = sum(tangent1 .* tangent2);  % from -1 to 1 with 0 being worst
    
%     obj = dist + derdiff;

obj = dist - dotprod^2;
end