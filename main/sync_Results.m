

local_results_dirs = {'C:\Users\rudi\Desktop\RD\Results\'  ,  'X:\Results\' , 'Y:\Results\' , 'Z:\Results\'};
% local_results_dirs = {'C:\Users\rudi\Desktop\RD\Results\'  , 'Z:\Results\'};
%  local_results_dirs = {'C:\Users\rudi\Desktop\RD\Results\' };

%  search_string = {'Results_efficiency_uber_refined_sweep1.mat' , 'Results_efficiency_uber_refined_sweep1_synced3.mat'};
% search_string = 'Results_Temporal_SN_radial_guesses_1.mat';
% search_string = {'Results_Dm_coarse_redo_1.mat' ,'Results_Dm_coarse_redo_2.mat', 'Results_Dm_refined_only_fixed.mat' , 'Results_Dm_few_redo_1.mat'};
% search_string = 'Results_Fore_Aft_SN.mat';
search_string = {'Results_Temporal_SN_radial_guesses_1_synced2.mat' , 'Results_Temporal_SN_cleanup_1.mat'};
    
search_string = {'Results_Temporal_SN_radial_guesses_1_synced2_cleaned.mat' , 'Results_Temporal_SN_cleanup_1.mat'};

search_string = {'Results_Temporal_SN_radial_guesses_1_synced2_cleanedagain.mat', 'Results_Temporal_SN_tryagain_1.mat'};

search_string = { 'Results_efficiency_pole2pole_*.mat'};
search_string = {'Results_Temporal_SN_shorter_max_arclength_1.mat'};

search_string = {'Results_Tumbling_Ease_SDE.mat'};

% search_string = {'Results_efficiency_refined_optimum_1_synced2_removed_1.45.mat', 'Results_efficiency_refined_optimum_1.mat'};

suffix = '_synced';  % append to end of new file names, use '' to overwrite originals

local_lock_namebase = 'lock';  %search for any files containing this

global_lock_file = 'C:\Users\rudi\Desktop\RD\Results\global_lock';

sync_freq = 60 * 90;  % seconds between syncs

timer = tic;
firstrun = true;

while true
    
    if toc(timer) < sync_freq && ~firstrun
        pause(10);
        continue
    end
    
    
    disp(['Attempting sync at ',datestr(now)]);
    
    clear fid
    worked = false;
    tic
    while ~worked
        
        
        
        
        locked = false;
        for i = 1:length(local_results_dirs)
            
            if ~isempty(dir([local_results_dirs{i}, '*',local_lock_namebase,'*']))
                %             if    exist( [local_results_dirs{i}, local_lock_filename],'file')
                locked = true;
                
                break
            end
            
        end
        
        if locked
            disp(['remote lock file found, paused for ',num2str(toc/60),' min']);
            pause(5);
        else
            
            worked = true;
            fid = fopen( [global_lock_file] , 'w' );  fclose(fid);
        end
        
    end
    
    global_Results = [];
    
    
    for dd = 1:length(local_results_dirs)
        
        if ~iscell(search_string)
            search_str = {search_string};
        else
            search_str = search_string;
        end
        temp = [];
        for ss = 1:length(search_str)
            temp = [temp;  dir([local_results_dirs{dd}, search_str{ss}])  ];  %could have multiple local sessions running at the same time with their own Results files
        end
        files = {temp.name};
        
        for f = 1:length(files)
            
            
            
            tic
            while true
                try
                    disp(['Loading ',[local_results_dirs{dd},files{f}] ]);
                    temp = load( [local_results_dirs{dd},files{f}]);
                    Results = temp.Results;
                    
                    empties = false(1,length(Results));
                    for rr = 1:length(Results)
                        if isempty(Results(rr).AR1)
                            empties(rr) = true;
                        end
                    end
                    Results(empties) = [];
                    
                    fields = fieldnames(Results);
                    break
                catch
                    disp(['can''t load remote results file, paused for ',num2str(toc/60),' min']);
                    pause(5)
                end
            end
            
            
            %%%%%%%%%%
            
            if isempty(global_Results)
                global_Results = Results;
                continue
            end
            
            global_geoms = [ [global_Results.AR1]'   [global_Results.AR2]'  [global_Results.amp]'  [global_Results.lambda]'  [global_Results.nlambda]'  ];
            local_geoms = [ [Results.AR1]'   [Results.AR2]'  [Results.amp]'  [Results.lambda]'  [Results.nlambda]'  ];
            
            
            for l = 1:length(Results)  %each geom result in local Results file
                
                
                
                [~,results_ind] =      ismember(local_geoms(l,:), global_geoms, 'rows');
                
                
                if results_ind == 0 % don't have this geom in Results yet, add a new entry
                    results_ind = length(global_Results) + 1;
                    % deal with local Results not having all the fields of
                    % global Results
                    names_global = fieldnames(global_Results);
                    names_local = fieldnames(Results);
                    
                    %                     if isequal( [local_results_dirs{dd},files{f}]  , 'Z:\Results\Results_1.mat') && l == 17834
                    %                         stopafra
                    %                     end
                    
                    if ~isequal(sort(names_global),sort(names_local))
                        temp = Results(l);
                        missing = setdiff(names_global, names_local);
                        for m = 1:length(missing)
                            temp.(missing{m}) = [];
                        end
                        global_Results(results_ind) = temp;
                    else % simple case, the fields match, just copy
                        if isequal((names_global),(names_local))
                            global_Results(results_ind) = Results(l);
                        else  %even though shatlab appeared to be cool with the two structures having the same fields in different order, it reliably crashes here unless we do it the long way:
                            for nl = 1:length(names_local)
                                global_Results(results_ind).(names_local{nl}) = Results(l).(names_local{nl});
                            end
                        end
                    end
                    
                    continue
                end
                
                fields = fieldnames(Results);
                for r = 1:length(fields)  %each field of this entry in local Results
                    
                    if ~isfield(global_Results(results_ind),fields{r})
                        global_Results(results_ind).(fields{r}) = Results(l).(fields{r});
                        continue
                    end
                    
                    if ~isempty(global_Results(results_ind).(fields{r})) && isempty(Results(l).(fields{r}))
                        continue  %don't overwrite actual results with empty placeholders
                    end
                    
                    global_Results(results_ind).(fields{r}) = Results(l).(fields{r});
                    
                end   %fields of each entry in local Results file
                
            end  % entries inside each local Results file
            
        end   % individual Results files inside each folder
        
    end  % local Results file folders
    
    
    
    for dd = 1:length(local_results_dirs)
        
        if ~iscell(search_string)
            search_str = {search_string};
        else
            search_str = search_string;
        end
        temp = [];
        for ss = 1:length(search_str)
            temp = [temp;  dir([local_results_dirs{dd}, search_str{ss}]) ];  %could have multiple local sessions running at the same time with their own Results files
        end
        files = {temp.name};
        
        for f = 1:length(files)
            
            tic
            while true
                try
                    Results = global_Results;
                    newname = [local_results_dirs{dd},files{f}(1:end-4),suffix,'.mat'];
                    disp(['Saving ',newname ]);
                    save( newname,'Results');
                    break
                catch
                    disp(['can''t save remote results file, paused for ',num2str(toc/60),' min']);
                    pause(5)
                end
            end
            break % only save one copy one time
        end
        break  % only save one copy one time
    end
    
    delete([global_lock_file]);
    
    disp(['Sync completed at ',datestr(now)]);
    
    timer = tic;  % restart countdown to next sync
    
    firstrun = false;
end %outer infinite loop that continuously runs






