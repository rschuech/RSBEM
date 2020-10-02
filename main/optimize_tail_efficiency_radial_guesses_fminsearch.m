% results_file = 'C:\Users\rudi\Desktop\RD\Results\results_opt_CFD01.mat';
rack = getenv('computername');

%   rack = 'CFD03';
% solutions_file = 'temporal_SN_solutions.mat'; %stores optimized tail parameters so far across all racks
% solutions_folders = {'C:\Users\rudi\Desktop\RD\Results\'  ,  'X:\Results\' , 'Y:\Results\' , 'Z:\Results\'}; % solutions file is only in one of them but path changes depending on rack
% solutions_folder = [];
% for s = 1:length(solutions_folders)
%     if exist( [solutions_folders{s},solutions_file] , 'file')
%         solutions_folder = solutions_folders{s};
%         break
%     end
% end

% if isempty(solutions_folder)
%     stopafra
% end
% others = setdiff({'CFD01','CFD02','CFD03', 'CFD04'},rack);
% for i = 1:length(others)
%     other_results_files{i} = ['C:\Users\rudi\Desktop\RD\Results\results_opt_',others{i},'.mat'];
% end

num = '1';

sweep_file = ['efficiency_smooth_CFD02.mat'];

clear X Fval Flag Output Time

temp = load(['C:\Users\rudi\Desktop\RD\sweeps\',sweep_file]);
bodies_sweep = temp.sweep.bodies;
guesses_sweep =[ temp.sweep.amp  temp.sweep.lambda temp.sweep.nlambda ];  %

scaling.mean = [0.4 5.5 1.42]';  %amp lambda nlambda
scaling.range = [0.4 4.6 0.1]';
scaling.obj = 3E4;
%%
% bodies_sweep = [10 0.4];
% guesses_sweep = [0.29 3.95 ];  % yields an arclength of almost 15 um
  
% max_arclength = 15;  %10 seems too small, an AR1 10 AR2 0 body has a definite effect, different than sphere or long curved body

%%

lock_file_name = ['results_lock_efficiency_',num];
results_file_name =   ['Results_efficiency_smooth_',num,'.mat'];
dump_folder = 'C:\Users\rudi\Desktop\RD\efficiency_dumps_smooth\';

% results_file_name = 'Results_.mat';

sweep_tempfile = ['tempsweep_efficiency_',num2str(num),'.mat'];

results_file = ['C:\Users\rudi\Desktop\RD\Results\',results_file_name];

%%
for sw = 1:114 %size(bodies_sweep,1)     %25:71
    
    AR1 = bodies_sweep(sw,1);  AR2 = bodies_sweep(sw,2);
    
    disp(['On opt sweep iter ',num2str(sw),' out of ',num2str(size(bodies_sweep,1)),'      body AR1 ',num2str(AR1),'  AR2 ',num2str(AR2)]);
    
    
     limits = [1 10; 0 1];
    guess = guesses_sweep(sw,:);
    
%     load([solutions_folder,solutions_file]); % solutions = [AR1 AR2 amp lambda nlambda]
%     if isempty(solutions) || sw == 1
%         guess  =   [0.29 3.95 ]; % halfway between sphere and AR1 10 AR2 0.4 optima
%     else
% %         solutions = solutions(1:end-1,:);
%         bodies_done = standardize(solutions(:,1:2),limits);
%         body = standardize([AR1 AR2],limits);
%         dists = sqrt(sum( (bodies_done - repmat(body,size(bodies_done,1),1)).^2 , 2) );
%         [~,ind] = min(dists);
%         guess = solutions(ind,3:4);
%     end
    

    
       lb = [1E-6   1E-6  0.1  ];
    ub = [6      20   3];
    
    if any(isnan(guess))
        guess = lb*3;  %prolly this is one of the donut bodies we couldn't guess a tail for
    end
    
%     if guess(1) < 0.3
%         guess = [0.4 3 1.3];
%     end
    
    objfun = @(tail) - 1/scaling.obj * optimization_wrapper_efficiency(AR1, AR2, tail(1)*scaling.range(1)+scaling.mean(1), tail(2)*scaling.range(2)+scaling.mean(2), tail(3)*scaling.range(3)+scaling.mean(3), results_file, lock_file_name, results_file_name , sweep_tempfile, dump_folder);
%     constraintfun = @(tail) max_arclength_constraint(tail,max_arclength);
    
    
%     optionsNM = optimset('display','iter','TolFun',1E6,'TolX',0.0005);
      optionsNM = optimset('display','iter','TolFun',1E6,'TolX',0.001);
    %     options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.01, 'AccelerateMesh',true,'maxfunctionevaluations',300, ...
    %         'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 0.5);
    % options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.25, 'AccelerateMesh',false,'maxfunctionevaluations',300, ...
    %         'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 0.5);
%     options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',2^(-7), 'AccelerateMesh',true,'maxfunctionevaluations',500, ...
%         'pollmethod','GPSPositiveBasis2N', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 1,'MeshTolerance',5E-4);
%       options = optimoptions('patternsearch','cache','on','cachetol',1E-6,'display','iter','InitialMeshSize',1, 'AccelerateMesh',true,'maxfunctionevaluations',500, ...
%         'pollmethod','GPSPositiveBasis2N',  'UseCompletePoll', false, 'MaxMeshSize' , 500,'MeshTolerance',0.0001,'PollOrderAlgorithm','Success','SearchFcn',{@searchneldermead,1,optionsNM});
    
%      options = optimoptions('patternsearch','cache','on','cachetol',1E-6,'display','iter','InitialMeshSize',0.01, 'AccelerateMesh',true,'maxfunctionevaluations',500, ...
%         'pollmethod','GPSPositiveBasis2N',  'UseCompletePoll', false, 'MaxMeshSize' , 500,'MeshTolerance',0.00002,'PollOrderAlgorithm','Success');
    
    
    
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
    
    output = objfun((guess(:) - scaling.mean)./scaling.range);
    
    continue
    % need to make sure intitial guess doesn't give NaN, which makes
    % patternsearch angry
    
    min_amp = lb(1);
    min_length = 0.03;  %min centerline length about same as tail radius, this would be a ridiculously short tail
    min_length = 0.01;
    
    while isnan(output)
        if guess(1) > min_amp %if we can reduce amp, try that
            disp('trying smaller amp')
            guess(1) = guess(1)*0.5
            output = objfun((guess(:) - scaling.mean)./scaling.range);
        elseif guess(3) * guess(2) > min_length  %if amp already tiny, try reducing tail centerline length too (via nlambda)
            disp('trying smaller nlambda')
            guess(3) = guess(3) * 0.5
            output = objfun((guess(:) - scaling.mean)./scaling.range);
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
    
    %PatternSearch
    %[X(:,sw), Fval(sw), Flag(sw), Output(sw)] = patternsearch(objfun,(guess(:) - scaling.mean)./scaling.range,[],[],[],[],(lb(:) - scaling.mean)./scaling.range,(ub(:) - scaling.mean)./scaling.range,[],options);

    %fminsearch
    [X(:,sw), Fval(sw), Flag(sw), Output(sw)] = fminsearch(objfun,(guess(:) - scaling.mean)./scaling.range,optionsNM);
%   

%[X(:,sw), Fval(sw), Flag(sw), Output(sw)] = fmincon(objfun,guess,[],[],[],[],lb,ub,[],options);
    % pause
    X(:,sw) = X(:,sw).*scaling.range + scaling.mean; 
    Fval(sw) = Fval(sw) * scaling.obj;
    X(:,sw)
    %Fval(sw)
    %     [X(:,sw), Fval(sw), Flag(sw), Output(sw)] = fmincon(objfun,guess,[],[],[],[],lb,ub,[],options)
    Time(sw) = toc(ptic) / 3600
    
    
%       load([solutions_folder,solutions_file]); % solutions = [AR1 AR2 amp lambda nlambda]
%     solutions(end+1,:) = [AR1 AR2 X(:,sw)'];
%     save([solutions_folder,solutions_file],'solutions');
    
    
end
%%

