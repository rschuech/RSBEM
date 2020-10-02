function out = crossprod(u,v)
% cross product with no error checking

out = [u(2)*v(3) - u(3)*v(2); u(3)*v(1) - u(1)*v(3); u(1)*v(2) - u(2)*v(1)];