% F_speed, F_curv are interpolants, takes standardized_AR1,
% standardized_AR2 as inputs

limits = [min(X(:)) max(X(:)); min(Y(:)) max(Y(:))];
[standardized_X, standardized_Y] = standardize(X,Y,limits);

obj_speed = @(AR) -F_speed(AR(1),AR(2));
obj_curv = @(AR) -F_curv(AR(1),AR(2));

speed_opt_pt = X_Y_unq(best_depvars == max(best_depvars),:);
[speed_opt_pt_standardized(1),speed_opt_pt_standardized(2)] = standardize(speed_opt_pt(1), speed_opt_pt(2), limits);
fmin_speed = obj_speed( speed_opt_pt_standardized );  % should be -1 for normalized power eff

curv_opt_pt = [X(1./mean_abs_curv == max(1./mean_abs_curv(:))) Y(1./mean_abs_curv == max(1./mean_abs_curv(:)))];
[curv_opt_pt_standardized(1),curv_opt_pt_standardized(2)] = standardize(curv_opt_pt(1), curv_opt_pt(2), limits);
fmin_curv = obj_curv( curv_opt_pt_standardized );

goal = [fmin_curv fmin_speed];

 fun = @(in) multiobj(in , obj_curv,  obj_speed);



nf = 2; % number of objective functions
N = 800; % number of points for plotting
onen = 1/N;
x = zeros(N+1,2);
f = zeros(N+1,nf);
x0 = [0.055 0.5];
options = optimoptions('fgoalattain','Display','off','MaxFunctionEvaluations',1000);
for r = 0:N
    r
    t = onen*r; % 0 through 1
    weight = [t,1-t];
    [x(r+1,:),f(r+1,:)] = fgoalattain(fun,x0,goal,weight,...
        [],[],[],[],[],[],[],options);
%     pause
end


%%

options = optimoptions('gamultiobj','PopulationSize',500,'FunctionTolerance',1e-6);
lb = [0 0];  ub = [1 0.8];
[x,f,exitflag] = gamultiobj(fun,2,[],[],[],[],lb,ub,options);

clear AR1 AR2
for i = 1:size(x,1)
[AR1(i), AR2(i)] = unstandardize(x(i,1),x(i,2),limits);
end


figure
plot(f(:,1),f(:,2),'k.');
xlabel('f_1')
ylabel('f_2')






standardized_ans = patternsearch(obj_speed,[0.2 0.5]',[],[],[],[],[0 0 ]',[1 1 ]' ,optimoptions('patternsearch','display','iter','meshtolerance',5E-8,'steptolerance',1E-10,'maxiterations',1200,'initialmeshsize',4))



