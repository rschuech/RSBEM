
    
    
loop_inds = randperm(length(phase_indices));
parfor(loop_ind = 1:length(phase_indices) , min([num_parallel_sessions, length(phase_indices)])   )
%     for loop_ind = 1:length(phase_indices) 

    warning('off','MATLAB:DELETE:FileNotFound');
    ti = phase_indices(loop_inds(loop_ind));
    
    ti_list_input = ti;
    
    pause( 30*rand);
    

   
    basename = ['hairs_',num2str(coplanar_length),'_',num2str(normal_length),'_step_',num2str(ti)];
    
        % if Parameters file exists, assume we already did this phase and skip it
    if ~isempty(   dir([folder_original , 'Parameters_step_',num2str(ti),'_*'])  )
        continue
    end
    
    %% modify template python scripts by inserting parameter values
    bat_file = [folder_temp,outputprefix,'_',basename,'_runSalome','.bat'];
    py_file = [folder_temp,outputprefix,'_',basename,'_runSalome','.py'];
    PID_file = [folder_temp,outputprefix,'_',basename,'_PID','.txt'];
%     kill_file = [folder_input,outputprefix,'_',basename,'_kill','.py'];
    
    %% bat file
    
    fid = fopen(run_bat_template);
    text = fread(fid, inf, '*char')';
    fclose(fid);
    
    text = regexprep(text, 'full_runSalome_path', py_file);
    
    if exist(bat_file,'file')
        delete(bat_file);
    end
    
    fid = fopen(bat_file, 'w');
    fwrite(fid, text);
    fclose(fid);
    %% run_Salome python file
    fid = fopen(run_py_template);
    text = fread(fid, inf, '*char')';
    fclose(fid);
    
    text = regexprep(text, 'PID_file', PID_file);
    
    if exist(py_file,'file')
        delete(py_file);
    end
    
    fid = fopen(py_file, 'w');
    fwrite(fid, text);
    fclose(fid);
    %% PID file is created by Salome/Python, just need to delete possible old file
    if exist(PID_file,'file')
        delete(PID_file);
    end
    
    %% main Salome python script
    
    
    fs = '%0.16g, ';  %formatting specification for writing values to text files
    
    
    template = [paths.infolder,paths.template_file];  
    
    fid = fopen(template);
    text = fread(fid, inf, '*char')';
    fclose(fid);
    
    
    python_file = [folder_temp,outputprefix,'_',basename,suffix,'.py'];
    if tail_mesh_only, temp = 'True'; else, temp = 'False'; end
    text = regexprep(text, 'tail_mesh_only_input', temp);
    text = regexprep(text, 'PID_file', PID_file);
    text = regexprep(text, 'folder_meshes', folder_original);
    text = regexprep(text, 'folder_temp', folder_temp);
    text = regexprep(text, 'parallel_hair_length_input', num2str(parallel_hair_length_input,digits));
    text = regexprep(text, 'normal_hair_length_top_input', num2str(normal_hair_length_top_input,digits));
    text = regexprep(text, 'normal_hair_length_bot_input', num2str(normal_hair_length_bot_input,digits));
    text = regexprep(text, 'normal_hairs_nsegs_radial_input', num2str(normal_hairs_nsegs_radial_input,digits));
    text = regexprep(text, 'parallel_hairs_nsegs_radial_input', num2str(parallel_hairs_nsegs_radial_input,digits));
    text = regexprep(text, 'ti_list_input', ['[',sprintf(fs,[ti_list_input]),']']);
    text = regexprep(text, '#_phase_pts', num2str(n_phase_pts));
    
    
    
    if exist(python_file,'file')
        delete(python_file);
    end
    
   
    
    fid = fopen(python_file, 'w');
    fwrite(fid, text);
    fclose(fid);
    
    
   
    %% do the deed - run Salome and create mesh (Salome saves the mesh .dat file itself)
    gar = [];  bage = [];
    switch getenv('computername')
        case 'UBERTOP'
            
            [gar,bage] = system([bat_file,' -t ',python_file]);
            
            
        case {'CFD01','CFD02','CFD03','CFD04'}
            
%             [gar,bage] = system(['C:\Users\rudi\Desktop\RD\SALOME-7.7.1-WIN64\run_salome.bat -t ',paths(1).python_file]);
            [gar,bage] = system([bat_file,' -t ',python_file]);
    end
    
   
    
    done_name = [folder_temp , 'done_step_' , num2str(ti) , '.txt'];
%     delete_files = {done_name, PID_file, kill_file,  [folder_input,outputprefix,'_',basename,'_runSalome','.*'] ,[folder_input,outputprefix,'_',basename,'.pyc']};
     delete_files = {done_name, PID_file,  [folder_input,outputprefix,'_',basename,'_runSalome','.*'] ,[folder_input,outputprefix,'_',basename,'.pyc']};
    
    times = 0:check_freq:max_time;
    
    for t = 1:length(times)
        
        pause(check_freq);
        
        if exist(done_name,'file')
            pause(25); %give Salome an extra long chance to close itself
            temp = load(PID_file);
            [~,~] = system(['taskkill /PID ',num2str(temp(1)), ' /T /F']); % the first PID is the main one
            for d = 1:length(delete_files)
                delete(delete_files{d});
            end
            break
        end
        
        if t == length(times)
            
           
            temp = load(PID_file);
            main_PID = temp(1);
            port = temp(end);
            
            
%             fid = fopen(kill_salome_carefully_template);
%             text = fread(fid, inf, '*char')';
%             fclose(fid);
            
%             text = regexprep(text, 'port_num', num2str(port));
            
%             if exist(kill_file,'file')
%                 delete(kill_file);
%             end
            
%             fid = fopen(kill_file, 'w');
%             fwrite(fid, text);
%             fclose(fid);
            
%               pause( 4 + (8 - 4)*rand);
%              [gar,bage] = system([bat_file,' -t ',kill_file]);
%              
%              
%              pause(25)
             
             
             [~,~] = system(['taskkill /PID ',num2str(main_PID),' /T /F']); % the first PID is the main one
             % by now, kill_file Salome should have recorded it's PID to the
             % PID_file
             temp = load(PID_file);
             kill_PID = temp(end-1);  % two new PIDs should have been appended by kill_file script, the first being the correct one like before
             [~,~] = system(['taskkill /PID ',num2str(kill_PID),' /T /F']); % the first PID is the main one
             
             for d = 1:length(delete_files)
                 delete(delete_files{d});
             end
        end
        
    end
    
    
end