

n = 7;
A = invhilb(n);
cond(A)

b = zeros(n,1);
%    b(3) = 1;
b(3) = 1;




dA = decomposition(A,'lu','CheckCondition',false);
sol = dA \ b;
sol


% attempt to reduce ill conditioning problem (if it is a problem) using
% iterative refinement with triple precision
% could do fancier job with symbolic toolbox ala Cleve Corner post but this
% takes 10 minutes to run.  Maybe other options using extended precision
% toolboxes or paid Advanpix toolbox
iterative_refinement = tic;
%     tol = 0.5;
%     improvement = Inf;
diff_norm_x = 1;  old_diff_norm_x = Inf;
while diff_norm_x < old_diff_norm_x && diff_norm_x > 1E-15
    r = residual3p(A,sol,b);
    % norm(r1) / (normA * norm(x))
    
    d = dA \ r;
    x_old = sol;
    sol = sol - d;
    old_diff_norm_x = diff_norm_x;
    diff_norm_x = norm(x_old - sol)
    
    
    sol
end