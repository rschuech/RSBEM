for i = 2:2:length(Inputs)  %only loop over timestepping runs
    i/length(Inputs)*2
    
    torque_file = [Inputs(i).paths.namebase.full,'_dump.mat'];
    timestepping_file = [Inputs(i).paths.namebase.body,'_',Inputs(i).paths.namebase.tail,'_motorBC_torque_timestepping.mat'];
    
    folder = Inputs(i).paths.dumpfolder;
    
    if xor(  exist([folder,torque_file]) , exist([folder,timestepping_file])  )  %found one but not other, problemo
        sstopafra
    elseif exist([folder,torque_file]) && exist([folder,timestepping_file])
        
        temp = load([folder,torque_file],'input');
        input = temp.input;
        
        Results(end+1).name = [input.paths.namebase.body, '_', input.paths.namebase.tail];
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