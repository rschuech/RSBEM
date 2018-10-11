function [u_vec] = compute_BCs_mex(forced_BCs,pts,refpoint,nthreads)

%applies translations and rotations defined in forced_BCs to the list of
%points in pts, typically a list of mesh vertices

%outputs u, a numel(pts) x 1 vector of velocity components at each pt in
%u v w order

%u = zeros(size(pts,1)*3,1);
velocity = NaN(size(pts,1),3);

parfor (i = 1:size(pts,1), nthreads)
    velocity(i,:) = forced_BCs.U + crossprod(forced_BCs.Omega, ( pts(i,:)' - refpoint ) );  %velocity of a pt = U + Omega cross R, R is vector from refpoint to current pt
end


velocity = velocity';  %switch to rows being x, y, z and cols being pts

u_vec = velocity(:);  %works since each column is placed under the one before, and each column is [u v w]' for a pt