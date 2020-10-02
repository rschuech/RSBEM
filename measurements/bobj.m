function [obj , k ] = bobj(inds,AR1AR2,shrink_factor)






subpts = AR1AR2(logical(inds),:);

[k,a] = boundary(subpts(:,1),subpts(:,2),shrink_factor);

[ geom ] = polygeom( subpts(k,1), subpts(k,2) ) ;  % k has last point same as first, which is OK for polygeom (whether included or not doesn't matter)
perim = geom(4);

obj = perim;