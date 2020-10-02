function [c,ceq] = bnonlcon(inds,AR1AR2,shrink_factor,min_count)

ceq = [];


subpts = AR1AR2(logical(inds),:);

[k] = boundary(subpts(:,1),subpts(:,2),shrink_factor);

in = inpolygon( AR1AR2(:,1),AR1AR2(:,2),subpts(k,1),subpts(k,2) );
% c = min_area - a;


c = min_count - sum(in);