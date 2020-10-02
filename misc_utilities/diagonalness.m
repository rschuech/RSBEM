function [d,tau_a] = diagonalness(angle, swimming_axis, D_aligned, varargin)

% if principle_axis_1_ind is included as optional input, it is used as primary principle axis index instead of the first column

d = NaN(size(angle)); tau_a = d;
   if isempty(varargin)
       inds = [2 3];
   else
 inds = setdiff(1:3,varargin{1});
   end
 
for i = 1:numel(angle)
    
    rotmat = rotate_arbitrary_vector(swimming_axis, angle(i));
    
    D_rot = rotmat * D_aligned * rotmat';
    
   d(i) = norm( D_rot - diag(diag(D_rot)) , 'fro' );
    
    
        tau_a(i) = 1/ (D_rot(inds(1),inds(1)) +  D_rot(inds(2),inds(2)));
        
end