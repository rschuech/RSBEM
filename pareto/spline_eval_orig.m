function [coord, tangent] = spline_eval(t, pp, ppder)

if isstruct(pp)
    coord = fnval(pp, t);
    temp = fnval(ppder, t);
    tangent = temp / sqrt(sum(temp.^2));
%     der = temp(2) / temp(1);
else
    coord = NaN;
%     der = NaN;
tangent = NaN;
end