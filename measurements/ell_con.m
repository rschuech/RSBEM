function [c,ceq,is_inside] = ell_con(ell, AR1AR2,min_count)


% ell = [a b x0 y0 angle]


ceq = [];

% https://stackoverflow.com/questions/7946187/point-and-ellipse-rotated-position-test-algorithm
is_inside = ( cos(ell(5)).*(AR1AR2(:,1) - ell(3)) + sin(ell(5)).*(AR1AR2(:,2) - ell(4)) ).^2  ./ ell(1)^2 + ...
    ( sin(ell(5)).*(AR1AR2(:,1) - ell(3)) - cos(ell(5)).*(AR1AR2(:,2) - ell(4)) ).^2  ./ ell(2).^2       <= 1;

c = min_count - sum(is_inside);


% subpts = AR1AR2(logical(inds),:);
% 
% [k] = boundary(subpts(:,1),subpts(:,2),shrink_factor);
% 
% in = inpolygon( AR1AR2(:,1),AR1AR2(:,2),subpts(k,1),subpts(k,2) );
% % c = min_area - a;
% 
% 
% c = min_count - sum(in);