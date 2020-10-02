function [obj, Curvature2] = objfun_curv(Curv, Curv2, constants, check_intersections, do_width_correction)

% Length_Curv = unstandardize( standardized_Length_Curv , limits);

% Length = Length_Curv(1);  Curvature = Length_Curv(2);

% Feret = Feret_Curv2(1);   Curv2 = Feret_Curv2(2);

[Feret, Curvature2, u_sol_best] = Length_Curv_2Feret(constants.Length, Curv, constants, check_intersections, do_width_correction);

% if Feret_Curv2(2) > 1E-2
%     obj = (  (  Feret - Feret_Curv2(1)  )  /  (  (Feret + Feret_Curv2(1)) / 2)   )^2   +     (  ( Curvature2 - Feret_Curv2(2) )  /  ((Curvature2 + Feret_Curv2(2)) / 2)   )^2;
% else
%     obj = (  (  Feret - Feret_Curv2(1)  )  /  (  (Feret + Feret_Curv2(1)) / 2)   )^2   +     (  ( Curvature2 - Feret_Curv2(2) )    )^2;
% end

% obj = (  (  Feret - Feret_Curv2(1)  )  /    Feret_Curv2(1)   )^2   +     (  ( Curvature2 - Feret_Curv2(2) ) / Feret_Curv2(2)   )^2;


% obj = abs(  (  Feret - Feret_Curv2(1)  )    )   +     abs(  ( Curvature2 - Feret_Curv2(2) )   );    % patternsearch

obj = Curvature2 - Curv2  ;

% function [standardized_coord] = standardize(coord, limits)

%coord is n X 2, n = number of points, col = AR1 AR2

% limits is [min_x max_x; min_y max_y], which can be beyond bounds of input
% x, y in case x, y is a subset of a larger parameter space


% standardized_Length_Curv = standardize( Length_Curv , limits);


