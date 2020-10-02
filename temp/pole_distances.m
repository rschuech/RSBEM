

clear constants

theta = Metadata.geom.nturns * 2 * pi;
constants.a = Metadata.geom.radius_curv * cos( theta - pi/2);
constants.b = Metadata.geom.radius_curv * sin( theta - pi/2);
constants.r = Metadata.geom.radius;
constants.x_1 = 0;  constants.y_1 = Metadata.geom.radius_curv;



guess = [4.5 2];
guess = [constants.x_1+constants.a*1.2   constants.y_1+constants.b*.9];
% guess = [-0.81 0.32];
% guess = [-0.75 0.36];

lb = [-Metadata.geom.radius_curv * 2  0];
ub = [Metadata.geom.radius_curv * 2   Metadata.geom.radius_curv * 2];

objfun = @(X) obj_fun_headpole(X,constants);

options = optimoptions('fmincon','StepTolerance',1E-16,'OptimalityTolerance',1E-12,'functiontolerance',1E-12);
answer = fmincon(objfun, guess, [],[],[],[],lb,ub,[],options)

% options = optimoptions('patternsearch','StepTolerance',1E-14,'functiontolerance',1E-12,'MeshTolerance',1E-13,'MaxIterations',1E6,'InitialMeshSize',0.5);
% answer = patternsearch(objfun, guess, [],[],[],[],lb,ub,[],options)

objfun(answer)

%   1.1614e-11