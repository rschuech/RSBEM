function [closest,ind] = near(x0, x)

x = x(:);

[~ , ind] = min(abs(x-x0));
closest = x(ind);