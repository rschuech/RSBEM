function [ p, idxs] = paretoFront( p )
% Filters a set of points P according to Pareto dominance, i.e., points
% that are dominated (both weakly and strongly) are filtered.
%
% Inputs: 
% - P    : N-by-D matrix, where N is the number of points and D is the 
%          number of elements (objectives) of each point.
%
% Outputs:
% - P    : Pareto-filtered P
% - idxs : indices of the non-dominated solutions
%
% Example:
% p = [1 1 1; 2 0 1; 2 -1 1; 1, 1, 0];
% [f, idxs] = paretoFront(p)
%     f = [1 1 1; 2 0 1]
%     idxs = [1; 2]

%  tol = 1E-10;  %how much better than equal does another point have to be to be "better"
% a_beats_b = @(a,b) (a - b)./((a+b)/2) >= tol;
% b_beats_a = @(a,b) (b - a)./((a+b)/2) >= tol;

 fun = @ge;
%  fun = a_beats_b;

[test_ind, ntasks] = size(p);  %test_ind is the pt currently being interrogated against all other remaining pts - starts at last point of input list
idxs = [1 : test_ind]';
while test_ind >= 1
    
    test_ind
    
    old_npts = size(p,1); % p will eventually just be optimal pts
%     indices = sum( bsxfun( @ge, p(i,:), p ), 2 ) == dim    &     any( bsxfun( fun, p(i,:), p ), 2 ) ;  
%   indices = any( bsxfun(a_beats_b, p(i,:), p ), 2 )    &     ~any( bsxfun( b_beats_a, p(i,:), p ), 2 ) ;  
    is_beaten = sum( bsxfun( fun, p(test_ind,:), p ), 2 ) == ntasks & any( bsxfun( @gt, p(test_ind,:), p ), 2) ;  % is current test_ind pt "better" at all tasks than each remaining pt?  is_beaten is length of remaining pts list, i.e. it refers to whether each remaining pt is beaten by the test point
    % originally ge, greater or equal, which is strong Pareto optimality
    % (beaten if another point is equal or better at each task)
    % the rare case of all tasks being equal is taken care of by 2nd check,
    % which makes sure the test point strictly beats the other point at at least one
    % task (i.e. if all tasks are equal, neither is removed)
%      indices = all( bsxfun(a_beats_b, p(i,:), p ), 2 ) ;
    is_beaten(test_ind) = false; % never count comparison of current test_ind pt to itself
   
    p(is_beaten,:) = []; %remove bad pts from current list
    idxs(is_beaten) = [];  % all points are optimal until proven suboptimal and removed 
    test_ind = test_ind - 1 - (old_npts - size(p,1)) + sum(is_beaten(test_ind:end));
end    

end