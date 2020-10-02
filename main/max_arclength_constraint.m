function [c,ceq] = max_arclength_constraint(tail,max_arclength)


c = bacterial_tail_arclength(tail) - max_arclength;  % should be <= 0

ceq = [];  

