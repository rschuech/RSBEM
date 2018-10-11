function [pole_separation] = line_projection_distance(pt1, pt2, slope)

% computes the distance between the projections of two points on a parametric line
% only the slope of the line is important (and only the direction, not
% magnitude of the slope vector)

% see Curved Rod Fore Aft Distance tablet note

 % doesn't matter whether head is in front of ass as usual or head is behind ass as in case of donuts


pole_separation = dot(  slope/sqrt(sum(slope.^2)) , pt2 - pt1 );


   