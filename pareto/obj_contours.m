function [obj] = obj_contours(t_0,  pp_0, F, limits)

obj = NaN(size(t_0));

for i = 1:numel(t_0)
    
    [coord_0] = spline_eval(t_0(i), pp_0);  % AR1 AR2 space
    [standardized_coord_0] = standardize(coord_0,limits);
    obj(i) = -F(standardized_coord_0(1),standardized_coord_0(2)); % attempt to maximize value of F interpolant along this contour (F is defined in standardized space)
end
% [obj] = objfun_contours(coord_0, coord_1, tangent_0, tangent_1);



