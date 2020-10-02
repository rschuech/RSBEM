function [Feret,u_sol_best] = Length2Feret(Length, constants, check_intersections, do_width_correction)

if do_width_correction
    constants.width = 2*Length/(pi - 4)*( sqrt( 1 + (pi - 4)/Length*constants.width ) - 1 );  % from mean width to actual width
end


constants.length = Length - constants.width;


[~,~,~,~,u1,u2,u3,shift] = curved_rod_parameters(constants);


if check_intersections
    u = linspace(0,1 - 1E-3,500);
    [x,y] = curved_rod_pts(u,constants,shift);
    if ~isempty(  intersections(x,y)  ) % self intersections detected
        Feret = NaN;
        return
    end
end
%%

u_poles = [u1/2  (u2+u3)/2];  % the poles
u_opp = u2;  % pt "opposite" u = 0 in the sense of the parameter u
% u_half = [(u1+u2)/2  (u3+1)/2]; % lower midpt, upper midpt
% u_guess = u_poles;
% u_guess = [0.4 .5];
lb = [0         u_opp];
ub = [u_opp     1    ];
n_guesses = 3;
u_guesses = [linspace(lb(1),ub(1),n_guesses)' linspace(lb(2),ub(2),n_guesses)'];

options = optimoptions('patternsearch', 'useparallel',false,'display','off');
clear u_sol feret
for g = 1:size(u_guesses,1)
    [u_sol(g,:),feret(g)] = patternsearch(@(X) obj_fun_feret(X,constants,shift), u_guesses(g,:), [],[],[],[],lb,ub,[],options);
end

[feret_best,ind] = min(feret);  feret_best = -feret_best;
u_sol_best = u_sol(ind,:);


Feret = feret_best;
