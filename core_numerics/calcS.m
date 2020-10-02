function [S] = calcS(x,x0,eps2)
%computes reg stokeslet S
%x is current location
%x0 is where reg stokeslet is
%eps2 is epsilon^2

d = x - x0;
%r = sqrt(sum(d.^2));
%r2 = sum(d.^2);
r2 = d(1)^2 + d(2)^2 + d(3)^2;
S = NaN(3,3);

c1 = (r2 + 2*eps2);
c2 =  (r2 + eps2)^(3/2);


% for i = 1:3
%     for j = 1:3
%         %
%         %       %  if i == j
% %         if i <= j
%             S(i,j) = ( (i == j)*c1 + d(i)*d(j) ) / c2;
% %         else
% %             S(i,j) = S(j,i);
% %         end
%         %       %  else
%         %          %   S(i,j) = d(i)*d(j) / (r^2 + epsilon^2)^(3/2);
%         %       %  end
%     end
% end

% small speed increase if everything is written out, but need to check mexed version

S(1,1) = ( c1 + d(1)^2 ) / c2;
S(1,2) = (      d(1)*d(2) ) / c2;
S(1,3) = (      d(1)*d(3) ) / c2;
S(2,1) = S(1,2);
S(2,2) = ( c1 + d(2)^2 ) / c2;
S(2,3) = (      d(2)*d(3) ) / c2;
S(3,1) = S(1,3);
S(3,2) = S(2,3);
S(3,3) = ( c1 + d(3)^2 ) / c2;


%S = (diag([c1 c1 c1]) + d*d') / ((r2 + eps2)^(3/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
