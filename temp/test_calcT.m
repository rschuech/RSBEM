function [] = test_calcT()

eps2 = 1E-2;

tic
for ii = 1:5E5
    x = rand(3,1);  x0 = rand(3,1);

d = x - x0;
%r = sqrt(sum(d.^2));
%r2 = sum(d.^2);
% r2 = d(1)^2 + d(2)^2 + d(3)^2;
T = NaN(3,3,3);

% c1 = (r2 + 2*eps2);
c1 =  (d(1)^2 + d(2)^2 + d(3)^2 + eps2)^(5/2);


for i = 1:3
    for j = 1:3
        for k = 1:3
%             if i <= j && j <= k
                T(i,j,k) = ( -6*d(i)*d(j)*d(k) - 3*eps2*( (j==k)*d(i) + (i==k)*d(j) + (i==j)*d(k) )  ) / c1;
%             else
%                 inds = sort([i j k]);
%                 T(i,j,k) = T(inds(1),inds(2),inds(3));
%             end
        end
    end
end


end


toc