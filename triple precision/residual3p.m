function r = residual3p(A,x,b)
% RESIDUAL3p Triple precision residual, A*x - b.
% r = residual3p(A,x,b) for matrix A, vectors x and b. 

   m = size(A,1);
   r = zeros(m,1);
   for k = 1:m
      r(k) = dot3p(A(k,:),x,-b(k));
   end

end
