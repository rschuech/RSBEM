function [mesh_succeed] = gen_mesh_wrapper(input)

max_mesh_time = 60 * 5;  %wait 5 minutes max before giving up waiting for Salome

factor = 1;  %initial tail mesh refinement factor from baseline parameters that usually work
min_factor = 0.7;

% factor = 0.45;
% min_factor = 0.1;

Amp = input.tail.amp;  Lambda = input.tail.lambda;  Nlambda = input.tail.nlambda;

disp(['Trying to generate tail with ','amp ',num2str(Amp,10),' lambda ',num2str(Lambda,10),' nlambda ',num2str(Nlambda,10),'    refinement factor ',num2str(factor)]);

gen_mesh;  %tries to mesh a tail with above parameters, also saves a corresponding metadata file

if isempty(geom) %we already created this mesh and it succeeded, done!
    if input.performance.verbose
    disp(['Already successfully generated tail with ','amp ',num2str(Amp,10),' lambda ',num2str(Lambda,10),' nlambda ',num2str(Nlambda,10)]);
    end
    mesh_succeed = true;
    return
end

pause(5);  %Salome is gonna take at least this long to do anything, prolly quite a lot longer


mesh_succeed = false;
if input.performance.verbose
disp('waiting for done file');
end
while ~ mesh_succeed
    
    
    
    mesh_tic = tic;
    salome_done = false;
    while toc(mesh_tic) < max_mesh_time
        
        if exist(paths.done_file,'file')  % it did a thing!
            if input.performance.verbose
            disp('done file found!');
            end
            delete(paths.done_file); %reset doneness signal
            salome_done = true;
            break
        else
            pause(2);  %give Salome this long to do something
        end
        
    end
    
    
    if salome_done
        if input.performance.verbose
        disp('checking mesh');
        end
        check_meshes;
        
        if succeeded
            mesh_succeed = true;
            break
        end
        
    end
    
    %at this point, either meshing clearly failed, or we never even saw a done
    %file.  Try refining and hope for the best....
    
    factor = factor - 0.05;  % increase refinement (decrease mesh sizes) by 5% each try
    
    if factor < min_factor
        disp('min refine factor reached; aborting')
        break
    end
    if input.performance.verbose
    disp(['Trying to generate tail with ','amp ',num2str(Amp,10),' lambda ',num2str(Lambda,10),' nlambda ',num2str(Nlambda,10),'    refinement factor ',num2str(factor)]);
    end
    gen_mesh;
    pause(5);
    
end

