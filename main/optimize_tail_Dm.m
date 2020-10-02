% results_file = 'C:\Users\rudi\Desktop\RD\Results\results_opt_CFD01.mat';
rack = getenv('computername');
results_file = ['C:\Users\rudi\Desktop\RD\Results\results_opt_Dm_',rack,'.mat'];

clear X Fval Flag Output Time

% temp = load(['C:\Users\rudi\Desktop\RD\sweep_',rack,'.mat']);
% bodies_sweep = temp.bodies;
% guesses_sweep = temp.guesses;


guesses_sweep = [0.2 3.5 1.3];
%%
for sw = 1:size(bodies_sweep,1)
    
    AR1 = bodies_sweep(sw,1);  AR2 = bodies_sweep(sw,2);
    
    disp(['On opt sweep iter ',num2str(sw),' out of ',num2str(size(bodies_sweep,1)),'      body AR1 ',num2str(AR1),'  AR2 ',num2str(AR2)]);
    
    
    
    guess = guesses_sweep(sw,:);
    
    lb = [0.2   3.5    1.3];
    ub = [1   3.5    1.3];
    
    
%         lb = [0.4   1    3];
%     ub = [0.4     8    3];
    
    objfun = @(tail) - optimization_wrapper_Dm(AR1, AR2, tail(1), tail(2), tail(3), results_file);
    
    
    
    
%     options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.01, 'AccelerateMesh',true,'maxfunctionevaluations',300, ...
%         'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 0.5);
% options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.25, 'AccelerateMesh',false,'maxfunctionevaluations',300, ...
%         'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 0.5);
options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.5, 'AccelerateMesh',true,'maxfunctionevaluations',500, ...
        'pollmethod','GPSPositiveBasis2N', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 0.5);
    
%      options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',1E-3, 'StepTolerance', 1E-3);
%       options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',1E-1, 'StepTolerance', 1E-3);
%        options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',2.5E-1, 'StepTolerance', 1E-2);
%        options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',1E-2, 'DiffMinChange', 1E-3, 'StepTolerance', 1E-3);
%               options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',2E-2, 'DiffMinChange', 5E-3, 'StepTolerance', 1E-3);
%            options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',5E-2, 'DiffMinChange', 1E-3, 'StepTolerance', 1E-3);
       
    ptic = tic;
    [X(:,sw), Fval(sw), Flag(sw), Output(sw)] = patternsearch(objfun,guess,[],[],[],[],lb,ub,[],options);
    X(:,sw)
    Fval(sw)
%     [X(:,sw), Fval(sw), Flag(sw), Output(sw)] = fmincon(objfun,guess,[],[],[],[],lb,ub,[],options)
    Time(sw) = toc(ptic) / 3600
    
    
   

    
end
%%

