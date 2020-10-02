Times = [];

cn = 0;
for nthreads = [1 2:2:20]
    cn = cn+1;
    
    Inputs.performance.nthreads = nthreads;
    close all
    timeall = tic;
    main;
    
    Times(cn) = toc(timeall);
    
end