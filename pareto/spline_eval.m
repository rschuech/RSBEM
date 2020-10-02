function [coord] = spline_eval(t, pp)

if isstruct(pp)
    coord = (fnval(pp, t))';

else
    coord = NaN;

end