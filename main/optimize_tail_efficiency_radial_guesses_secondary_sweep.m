% results_file = 'C:\Users\rudi\Desktop\RD\Results\results_opt_CFD01.mat';
rack = getenv('computername');

%   rack = 'CFD03';
% solutions_file = 'solutions.mat'; %stores optimized tail parameters so far across all racks
% solutions_folders = {'C:\Users\rudi\Desktop\RD\Results\'  ,  'X:\Results\' , 'Y:\Results\' , 'Z:\Results\'}; % solutions file is only in one of them but path changes depending on rack
% solutions_folder = [];
% for s = 1:length(solutions_folders)
%     if exist( [solutions_folders{s},solutions_file] , 'file')
%         solutions_folder = solutions_folders{s};
%         break
%     end
% end
% 
% if isempty(solutions_folder)
%     stopafra
% end

% others = setdiff({'CFD01','CFD02','CFD03', 'CFD04'},rack);
% for i = 1:length(others)
%     other_results_files{i} = ['C:\Users\rudi\Desktop\RD\Results\results_opt_',others{i},'.mat'];
% end

num = '1';

sweep_file = ['secondary_sweep_efficiency_radially_ordered_',rack,'_',num2str(num),'.mat'];

clear X Fval Flag Output Time

temp = load(['C:\Users\rudi\Desktop\RD\sweeps\',sweep_file]);
sweep_subset = temp.sweep_subset;
bodies_sweep = temp.sweep_subset(:,1:2);
% guesses_sweep = temp.guesses_combined;
dump_folder = 'C:\Users\rudi\Desktop\RD\efficiency_dumps_radial_guesses\';

lock_file_name = ['results_lock_',num];
results_file_name =   ['Results_efficiency_secondary_sweep',num,'.mat'];

% results_file_name = 'Results_.mat';

sweep_tempfile = ['tempsweep_',num2str(num),'.mat'];

results_file = ['C:\Users\rudi\Desktop\RD\Results\',results_file_name];
%CFD03   45:71    82:end
%CFD02  56:end  CFD03
% CFD04  61:77    88:end   35:43
% CFD01 35:47
%%
for sw = 1:size(bodies_sweep,1)     %25:71
    
    AR1 = bodies_sweep(sw,1);  AR2 = bodies_sweep(sw,2);
    
    disp(['On opt sweep iter ',num2str(sw),' out of ',num2str(size(bodies_sweep,1)),'      body AR1 ',num2str(AR1),'  AR2 ',num2str(AR2)]);
    
    
%     limits = [1 10; 0 1];
% %     guess = guesses_sweep(sw,:);
%     
%     load([solutions_folder,solutions_file]); % solutions = [AR1 AR2 amp lambda nlambda]
%     if isempty(solutions)
%         guess  =    [0.6    4.2     1.35   ];
%     else
%         bodies_done = standardize(solutions(:,1:2),limits);
%         body = standardize([AR1 AR2],limits);
%         dists = sqrt(sum( (bodies_done - repmat(body,size(bodies_done,1),1)).^2 , 2) );
%         [~,ind] = min(dists);
%         guess = solutions(ind,3:5);
%     end

guess = sweep_subset(sw,3:5);
        
    lb = [0.001   0.1    0.05];
    ub = [6      10     8];
    
    if any(isnan(guess))
        guess = lb*3;  %prolly this is one of the donut bodies we couldn't guess a tail for
    end
    
%     if guess(1) < 0.3
%         guess = [0.4 3 1.3];
%     end
    
    objfun = @(tail) - optimization_wrapper_redo_timestepping_medium(AR1, AR2, tail(1), tail(2), tail(3), results_file, lock_file_name, results_file_name , sweep_tempfile, dump_folder);
    
    
    
    
    %     options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.01, 'AccelerateMesh',true,'maxfunctionevaluations',300, ...
    %         'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 0.5);
    % options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',0.25, 'AccelerateMesh',false,'maxfunctionevaluations',300, ...
    %         'pollmethod','GPSPositiveBasisNp1', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 0.5);
%     options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',2^(-7), 'AccelerateMesh',true,'maxfunctionevaluations',500, ...
%         'pollmethod','GPSPositiveBasis2N', 'PollOrderAlgorithm', 'random', 'UseCompletePoll', false, 'MaxMeshSize' , 1,'MeshTolerance',5E-4);
      options = optimoptions('patternsearch','cache','on','cachetol',1E-4,'display','iter','InitialMeshSize',2^(-3), 'AccelerateMesh',true,'maxfunctionevaluations',500, ...
        'pollmethod','GPSPositiveBasis2N',  'UseCompletePoll', false, 'MaxMeshSize' , 8,'MeshTolerance',0.005,'PollOrderAlgorithm','Success');
    
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
    
    
%     load([solutions_folder,solutions_file]); % solutions = [AR1 AR2 amp lambda nlambda]
%     solutions(end+1,:) = [AR1 AR2 X(:,sw)'];
%     save([solutions_folder,solutions_file],'solutions');
    
    
end
%%

