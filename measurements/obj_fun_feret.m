function [obj] = obj_fun_feret(u,constants,shift)

distance = @(x,y) sqrt( diff(x,[],2).^2 + diff(y,[],2).^2 );

[x,y] = curved_rod_pts(u,constants,shift);

obj = -distance( x,y );