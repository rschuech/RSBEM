function [T] = calcT(x,x0,eps2)
%computes reg stresslet T for double layer integrals
%x is current location
%x0 is where reg stresslet is
%eps2 is epsilon^2

d = x - x0;
%r = sqrt(sum(d.^2));
%r2 = sum(d.^2);
% r2 = d(1)^2 + d(2)^2 + d(3)^2;
T = NaN(3,3,3);

% c1 = (r2 + 2*eps2);
c1 =  (d(1)^2 + d(2)^2 + d(3)^2 + eps2)^(5/2);

% somehow, without shortcut is 5 times faster unmexed, but mexed might be
% different - need to check

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

% writing everything out like calcS may be fastest?  need to check mexed
% versions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
