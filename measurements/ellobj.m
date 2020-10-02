function [obj]  = ellobj(ell)


% ell = [a b x0 y0 angle]

obj = pi*ell(1)*ell(2);

% subpts = AR1AR2(logical(inds),:);
% 
% [k,a] = boundary(subpts(:,1),subpts(:,2),shrink_factor);
% 
% [ geom ] = polygeom( subpts(k,1), subpts(k,2) ) ;  % k has last point same as first, which is OK for polygeom (whether included or not doesn't matter)
% perim = geom(4);
% 
% obj = perim;