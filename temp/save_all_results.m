folder = '../swept_dumps/';


files = dir([folder,'*torque_dump*.mat']);
files = {files.name};




for i = 1:length(files)  %only loop over timestepping runs
    i/length(files)
    
    torque_file = files{i};
    timestepping_file = [files{i}(1:end-8),'timestepping.mat'];
    

    
    if xor(  exist([folder,torque_file]) , exist([folder,timestepping_file])  )  %found one but not other, problemo
        sstopafra
    elseif exist([folder,torque_file]) && exist([folder,timestepping_file])
        
        temp = load([folder,torque_file],'input');
        input = temp.input;
        try
        Results(end+1).name = [input.paths.namebase.body, '_', input.paths.namebase.tail];
        catch
             Results(end+1).name = [input.paths.namebase_body, '_', input.paths.namebase_tail];
        end
        Results(end).AR1 = input.body.AR(1);
        Results(end).AR2 = input.body.AR(2);
        Results(end).amp = input.tail.amp;
        Results(end).lambda = input.tail.lambda;
        Results(end).nlambda = input.tail.nlambda;
        
        temp = load([folder,timestepping_file],'fits','avg_omega','motor_torque');
        
        Results(end).Avg_Speed = temp.fits.converged.speed;
        Results(end).diff_cutoff = temp.fits.converged.diff_cutoff;
        Results(end).t_convergence = temp.fits.T( temp.fits.line.speed == temp.fits.converged.speed);
        Results(end).Avg_Omega = temp.avg_omega;
        Results(end).Motor_Torque = temp.motor_torque;
    end
    
end
%                 save(input.paths.results_file,'Results');