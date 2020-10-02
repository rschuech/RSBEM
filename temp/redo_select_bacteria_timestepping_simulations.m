

folder = 'C:\Users\rudi\Desktop\RD\swimming dumps\';
old_folder = 'C:\Users\rudi\Desktop\RD\swimming dumps\old broken\';
interp_dumps = dir([folder,'*dump.mat']);
interp_dumps = {interp_dumps.name};

for dumpi = 1:length(interp_dumps)
    dumpi
    load([folder,interp_dumps{dumpi}]);
    temp = load(  [old_folder,input.paths.namebase.full,'_timestepping','.mat']  );
    
    input.paths.dumpfolder = 'C:\Users\rudi\Desktop\RD\swimming dumps\';
    input.performance.verbose = true;
    
    dump.y0 = y0;  dump.best_interpolant = best_interpolant;
   
    
    timestepping_solution0 = [];  refpoint0 = Mesh(1).refpoints(:,1);
    [fits, timestepping_solution, timings] = timestepping(input,dump,timestepping_solution0,refpoint0);  %timestepping dump file is saved internally, contains timestepping_solution and fits
    %             if input.performance.verbose
    
    interpolant = interpolant(end).omega;  %last and most accurate interpolant for tail rotation rate omega
    
    [avg_omega] = compute_avg_omega(interpolant);
    
    %in addition to storing in memory for later inclusion into aggregate results file, immediately save avg_omega and motor_torque inside timestepping dump file
                m = matfile([input.paths.dumpfolder,input.paths.namebase.full,'_timestepping','.mat'],'Writable',true);
                m.avg_omega = avg_omega;
                m.motor_torque = input.tail.motor_torque;
    
                
               [ fits.converged.speed    temp.fits.converged.speed   (fits.converged.speed -   temp.fits.converged.speed ) / fits.converged.speed * 100]
                
end


