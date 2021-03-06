% results_file = 'C:\Users\rudi\Desktop\RD\Results\results_opt_CFD01.mat';
rack = getenv('computername');

%  rack4 = 'CFD02';

% others = setdiff({'CFD01','CFD02','CFD03', 'CFD04'},rack);
% for i = 1:length(others)
%     other_results_files{i} = ['C:\Users\rudi\Desktop\RD\Results\results_opt_',others{i},'.mat'];
% end

num = '1';

sweep_file = ['sweep_fore_aft_SN_',rack,'_',num2str(num),'.mat'];

clear X Fval Flag Output Time

% temp = load(['C:\Users\rudi\Desktop\RD\sweeps\',sweep_file]);
% bodies_sweep = temp.bodies;
% guesses_sweep = temp.guesses;
%%
bodies_sweep = [5 0.5];
% guesses_sweep = [0.1 3 ];  % yields an arclength of almost 15 um

max_arclength = 15;  %10 seems too small, an AR1 10 AR2 0 body has a definite effect, different than sphere or long curved body

%%

lock_file_name = ['results_lock_fore_aft_SN_',num];
results_file_name =   ['Results_fore_aft_SN_',num,'.mat'];

% results_file_name = 'Results_.mat';

sweep_tempfile = ['tempsweep_fore_aft_SN_',num2str(num),'.mat'];

results_file = ['C:\Users\rudi\Desktop\RD\Results\',results_file_name];
%CFD03   45:71    82:end
%CFD02  56:end  CFD03
% CFD04  61:77    88:end   35:43
% CFD01 35:47
%%
for sw = 1:size(bodies_sweep,1)     %25:71
    
    AR1 = bodies_sweep(sw,1);  AR2 = bodies_sweep(sw,2);
    
    disp(['On opt sweep iter ',num2str(sw),' out of ',num2str(size(bodies_sweep,1)),'      body AR1 ',num2str(AR1),'  AR2 ',num2str(AR2)]);
    
    
    
    %     guess = guesses_sweep(sw,:);
    guess = [0.1 5 ]; % halfway between sphere and AR1 10 AR2 0.4 optima
    %     lb = [1E-6   1E-6    1E-6];
    %     ub = [6      20     10];
    
    lb = [0  0    ];
    ub = [6      5   ];
    
    if any(isnan(guess))
        guess = lb*3;  %prolly this is one of the donut bodies we couldn't guess a tail for
    end
    
    %     if guess(1) < 0.3
    %         guess = [0.4 3 1.3];
    %     end
    
    objfun = @(tail) - optimization_wrapper_fore_aft_SN(AR1, AR2, tail(1), tail(2), max_arclength, results_file, lock_file_name, results_file_name , sweep_tempfile);
    %     constraintfun = @(tail) max_arclength_constraint(tail,max_arclength);
    
    
    
    %     options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.01, 'AccelerateMesh',true,'maxfunctionevaluations',300, ...
    %         'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 0.5);
    % options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.25, 'AccelerateMesh',false,'maxfunctionevaluations',300, ...
    %         'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 0.5);
    %     options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',2^(-7), 'AccelerateMesh',true,'maxfunctionevaluations',500, ...
    %         'pollmethod','GPSPositiveBasis2N', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 1,'MeshTolerance',5E-4);
    options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',20, 'AccelerateMesh',true,'maxfunctionevaluations',500, ...
        'pollmethod','GPSPositiveBasis2N',  'UseCompletePoll', false, 'MaxMeshSize' , 500,'MeshTolerance',0.5,'PollOrderAlgorithm','Success');
    
    %      options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',1E-3, 'StepTolerance', 1E-3);
    %       options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',1E-1, 'StepTolerance', 1E-3);
    %        options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',2.5E-1, 'StepTolerance', 1E-2);
    %        options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',1E-2, 'DiffMinChange', 1E-3, 'StepTolerance', 1E-3);
    %               options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',2E-2, 'DiffMinChange', 5E-3, 'StepTolerance', 1E-3);
    %            options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',2E-5,'FiniteDifferenceStepSize',5E-2, 'DiffMinChange', 1E-3, 'StepTolerance', 1E-3);
    %     options = optimoptions('fmincon','display','iter-detailed','FunctionTolerance',1E-4,'FiniteDifferenceStepSize',1E-3, 'DiffMinChange', 1E-5, 'StepTolerance', 2E-4);
    % FunctionTolerance not actually used by default interior-point
    % algorithm
    
    ptic = tic;
    
    output = objfun(guess);
    % need to make sure intitial guess doesn't give NaN, which makes
    % patternsearch angry
    
    min_amp = lb(1);
    min_length = 0.03;  %min centerline length about same as tail radius, this would be a ridiculously short tail
    min_length = 0.01;
    
    while isnan(output)
        if guess(1) > min_amp %if we can reduce amp, try that
            disp('trying smaller amp')
            guess(1) = guess(1)*0.5
            output = objfun(guess);
        elseif guess(3) * guess(2) > min_length  %if amp already tiny, try reducing tail centerline length too (via nlambda)
            disp('trying smaller nlambda')
            guess(3) = guess(3) * 0.5
            output = objfun(guess);
        else %give up
            disp('giving up')
            output = NaN;
            break;
        end
        
    end
    
    %at this point, guess should yield non-NaN output
    
    if isnan(output)
        %X(:,sw) = NaN; Fval(sw) = NaN; Flag(sw) = NaN; Output(sw) = NaN;
        continue
    end
    
    [X(:,sw), Fval(sw), Flag(sw), Output(sw)] = patternsearch(objfun,guess,[],[],[],[],lb,ub,[],options);
    %     [X(:,sw), Fval(sw), Flag(sw), Output(sw)] = fmincon(objfun,guess,[],[],[],[],lb,ub,[],options);
    % pause
    
    X(:,sw)
    Fval(sw)
    %     [X(:,sw), Fval(sw), Flag(sw), Output(sw)] = fmincon(objfun,guess,[],[],[],[],lb,ub,[],options)
    Time(sw) = toc(ptic) / 3600
    
    
    
    
    
end
%%

%% test of different amps for arclength = 15, lambda = 5
% indeed monotonically increasing tau_a with decreasing amp, right down to
% 0 amp
% clear SN
tails = [0.1 5; 0.05 5; 0.025 5; 0.01 5;   0.005 5; 0.0025 5; 0.001 5;  0 5;];
% tails = [0.005 5; 0.0025 5; 0.001 5;];
for sw = 8:size(tails,1)
    sw
    SN(sw) = objfun(tails(sw,:)');
end
