function [coord] = unstandardize(standardized_coord,limits)

%standardized_coord is n X 2, n = number of points, col = AR1 AR2

% limits is [min_x max_x; min_y max_y], which can be beyond bounds of input
% x, y in case x, y is a subset of a larger parameter space

coord(:,1) = (limits(1,2) - limits(1,1)) * standardized_coord(:,1) + limits(1,1);
coord(:,2) = (limits(2,2) - limits(2,1)) * standardized_coord(:,2) + limits(2,1);