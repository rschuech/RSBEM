% Inputs = Inputs(1:14);
% Inputs = Inputs(1:196);
%  Inputs = Inputs(31:46);
% Inputs = Inputs(31:60);
% Inputs = Inputs(61:74);
% Inputs = Inputs(61:90);
% Inputs = Inputs(197:394);
Inputs = Inputs(395:592);

Inputs = Inputs(randperm(length(Inputs)));
% do freeswimming timestepping afterwards in parallel to maximize parallel scaling,
% since each timestepping run can only use 4 CPUs max when run one at a time

problemtypes = {Inputs.problemtype};
is_freeswim = strcmp(problemtypes,'freeswim');
%inds = find(is_freeswim & ~bad_sweep_i);  %must be a freeswim case, and can't have been a skipped run due to no mesh file
inds = find(is_freeswim );  %must be a freeswim case, and can't have been a skipped run due to no mesh file

if isempty(inds)
    return
end

init_parpool(Inputs(1).performance.nthreads, length(inds));


parfor (sweep_i = 1:length(inds))
    
    %for sweep_i = 1:length(inds)
    ind = inds(sweep_i);
    input = Inputs(ind);
    if exist([input.paths.dumpfolder,input.paths.namebase.full,'_dump','.mat'],'file') && ~exist([input.paths.dumpfolder,input.paths.namebase.full,'_timestepping','.mat'],'file')
              dump = load([input.paths.dumpfolder,input.paths.namebase.full,'_dump','.mat'],'y0','best_interpolant');
%         try
            timestepping(input,dump);
%         catch
%             disp(['Problem with iter ',num2str(sweep_i)]);
%             pause
%         end
        
    end
    
    
    
    % disp([input.paths.fullnamebase,'    timestepping took ',num2str(toc(timestepping_tic)/60), ' min']);
    
    
end