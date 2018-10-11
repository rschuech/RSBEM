function [c,ceq] = unit_slope_constraint(slope)

c = [];

ceq = sqrt( sum( slope.^2) ) - 1;