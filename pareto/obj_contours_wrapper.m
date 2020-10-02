function [obj] = obj_contours_wrapper(t_0, c_1, t_1,  pp_0, ppder_0, x_1, y_1, Z_1)

[ contours_1] = get_contours(c_1, x_1,y_1,Z_1);


[pp_1, ppder_1 ] = splining(contours_1);


[coord_0, tangent_0] = spline_eval(t_0, pp_0, ppder_0);
[coord_1, tangent_1] = spline_eval(t_1, pp_1, ppder_1);

[obj] = objfun_contours(coord_0, coord_1, tangent_0, tangent_1);



