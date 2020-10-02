function [Mesh] = compute_BCs_mex(resistance_BCs,Mesh,assembly_input)

%calculates velocities of each vert in each submesh of Mesh resulting from
% rigid body translational and rotational velocities defined in each element of resistance_BCs to

%adds field u to Mesh, velocities at each vertex separate for each submesh

[refpoint] = get_rotational_references(Mesh, assembly_input);

for f = 1:length(resistance_BCs)
    
    
    for m = 1:length(Mesh)
        
        % 3rd dim is for each BC, e.g. x y z translation, x y z rotation, 6 total for diffusion tensor calcs but prolly just 1 BC for mobility problems
        
        vel = NaN(Mesh(m).n_nodes,3); % need temp variable to avoid parfor error
        parfor (i = 1:Mesh(m).n_nodes, assembly_input.performance.nthreads)
            vel(i,:) = resistance_BCs(f).U + crossprod(resistance_BCs(f).Omega, ( Mesh(m).nodes(i,:)' - refpoint ) );  %velocity of a pt = U + Omega cross R, R is vector from refpoint to current pt
        end
        Mesh(m).u(:,:,f) = vel;
        
    end
    
end