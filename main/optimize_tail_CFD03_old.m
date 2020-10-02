% results_file = 'C:\Users\rudi\Desktop\RD\Results\results_opt_CFD01.mat';
 rack = getenv('computername');
results_file = ['C:\Users\rudi\Desktop\RD\Results\results_opt_',rack,'.mat'];

% m = matfile('C:\Users\rudi\Desktop\RD\Results\patternsearch_test.mat','writable',true);
% m.X{1:3} = X;

AR1 = 2.5;  AR2 = 0.25;

% guess = [0.43 3 1.4]';

guess = [0.45643125000000  3.20613125000000  1.309503125000000]'  + 0.001;


guess = [ 0.458837500000000    3.232131250000000     1.314878125000000  ]';

lb = guess - guess * 0.1;

ub = guess + guess*0.1;

objfun = @(tail) - optimization_wrapper(AR1, AR2, tail(1), tail(2), tail(3), results_file);
%% fmincon

options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',1E-3, 'StepTolerance', 1E-3);
[x, fval, flag, output] = fmincon(objfun,guess,[],[],[],[],lb,ub,[],options)
% 

%%
  options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.01, 'AccelerateMesh',true,'maxfunctionevaluations',300, ...
        'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false);
    
    [x, fval, flag, output] = patternsearch(objfun,guess,[],[],[],[],lb,ub,[],options)
    
%%

%%
clear options

options{1} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasis2N', 'PollOrderAlgorithm', 'consecutive', 'UseCompletePoll', false);

options{2} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasis2N', 'PollOrderAlgorithm', 'success', 'UseCompletePoll', false);

options{3} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasis2N', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false);

options{4} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasis2N', 'PollOrderAlgorithm', 'consecutive', 'UseCompletePoll', true);

options{5} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasis2N', 'PollOrderAlgorithm', 'success', 'UseCompletePoll', true);

options{6} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasis2N', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', true);

%%
options{7} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'consecutive', 'UseCompletePoll', false);

options{8} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'success', 'UseCompletePoll', false);

options{9} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.0001, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false);

options{10} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'consecutive', 'UseCompletePoll', true);

options{11} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'success', 'UseCompletePoll', true);

options{12} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', true);
%%
options{13} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasis2N', 'PollOrderAlgorithm', 'consecutive', 'UseCompletePoll', false);

options{14} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasis2N', 'PollOrderAlgorithm', 'success', 'UseCompletePoll', false);

options{15} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasis2N', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false);

options{16} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasis2N', 'PollOrderAlgorithm', 'consecutive', 'UseCompletePoll', true);

options{17} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasis2N', 'PollOrderAlgorithm', 'success', 'UseCompletePoll', true);

options{18} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasis2N', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', true);
%%
options{19} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasisNp1', 'PollOrderAlgorithm', 'consecutive', 'UseCompletePoll', false);

options{20} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasisNp1', 'PollOrderAlgorithm', 'success', 'UseCompletePoll', false);

options{21} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false);

options{22} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasisNp1', 'PollOrderAlgorithm', 'consecutive', 'UseCompletePoll', true);

options{23} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasisNp1', 'PollOrderAlgorithm', 'success', 'UseCompletePoll', true);

options{24} = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.005, 'AccelerateMesh',true,'maxfunctionevaluations',150, ...
    'pollmethod','MADSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', true);

clear X Fval Flag Output
for i = 9

ptic = tic;
[x, fval, flag, output] = patternsearch(objfun,guess,[],[],[],[],lb,ub,[],options{i})
toc(ptic)

X{i} = x;
Fval{i} = fval;
Flag{i} = flag;
Output{i} = output;

end

% fval =
% 
%   -0.006383079895653
% 
% x
% 
% x =
% 
%    0.461181250000000
%    3.247131250000000
%    1.313003125000000
%%

%    34      136    -0.00636551      0.003906     Refine Mesh
%    35      136    -0.00636551      0.001953     Refine Mesh
%    36      136    -0.00636551     0.0009766     Refine Mesh
%    37      136    -0.00636551     0.0004883     Refine Mesh
%    38      136    -0.00636551     0.0002441     Refine Mesh
%    39      136    -0.00636551     0.0001221     Refine Mesh
%    40      136    -0.00636551     6.104e-05     Refine Mesh
%    41      136    -0.00636551     1.526e-05     Refine Mesh
%    42      136    -0.00636551     1.907e-06     Refine Mesh
%    43      136    -0.00636551     1.192e-07     Refine Mesh
% Optimization terminated: mesh size less than options.MeshTolerance.
% 
% x =
% 
%    0.456431250000000
%    3.206131250000000
%    1.309503125000000
% 
% 
% fval =
% 
%   -0.006365511696556
% 
% 
% flag =
% 
%      1
% 
% 
% output = 
% 
%          function: @(tail)-optimization_wrapper(2,0.2,tail(1),tail(2),tail(3),results_file)
%       problemtype: 'boundconstraints'
%        pollmethod: 'gpspositivebasis2n'
%     maxconstraint: 0
%      searchmethod: []
%        iterations: 43
%         funccount: 136
%          meshsize: 1.192092895507813e-07
%          rngstate: [1x1 struct]
%           message: 'Optimization terminated: mesh size less than options.MeshTolerance.'