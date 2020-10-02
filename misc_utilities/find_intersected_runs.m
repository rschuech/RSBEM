 
% folder = 'C:\Users\rudi\Desktop\RD\swept_dumps_archived\';
folder = 'C:\Users\rudi\Desktop\RD\swept_dumps\';


intersected_files = {};

files = dir([folder,'*torque_dump.mat']);
files = {files.name};



ppm = ParforProgressStarter2('shat', length(files));

parfor f = 1:length(files)  %timestepping dumps
    dump_file = files{f};
    
    
    temp = load([folder,dump_file]);
    
    input = temp.input;
    Mesh = temp.Mesh;
    
    
    if isfield(input.paths,'namebase')
        name = input.paths.namebase.full;
    else
        name = input.paths.fullnamebase;
    end
    
    
    
    
    if isfield(input.accuracy, 'check_intersections_tolerance') && input.accuracy.check_intersections_tolerance == input.tail.radius / 2 && input.accuracy.check_intersections_n_angles >= 100
        ppm.increment(f);
        continue
    end
    
    
    
    
    input.accuracy.check_intersections_tolerance = input.tail.radius / 2;
    input.accuracy.check_intersections_n_angles = 100;
    
    

    [is_intersected] = submesh_intersections(Mesh,input.accuracy.check_intersections_tolerance, input.accuracy.check_intersections_n_angles,false,1);  % checks for self-intersection at each angle in parallel

    
    if is_intersected
        
        if isfield(input.paths,'namebase')
            
            namebase = [input.paths.namebase.body,'_',input.paths.namebase.tail,'_'];
        else
            namebase = [input.paths.namebase_body,'_',input.paths.namebase_tail,'_'];
        end
        
        bad_files = dir([folder,namebase,'*']);
        bad_files = {bad_files.name};
        
        intersected_files = [intersected_files; bad_files(:)];
        
    end
    
    ppm.increment(f);
    
    
    
end

delete(ppm);
