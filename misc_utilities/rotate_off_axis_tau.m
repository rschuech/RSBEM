function [tau_a] = rotate_off_axis_tau(angle, swimming_axis, D_aligned, varargin)

% if principle_axis_1_ind is included as optional input, it is used as primary principle axis index instead of the first column

tau_a = NaN(size(angle));

for i = 1:numel(angle)
    
    rotmat = rotate_arbitrary_vector(swimming_axis, angle(i));
    
    D_rot = rotmat * D_aligned * rotmat';
    
    if isempty(varargin)
        tau_a(i) = 1/ (D_rot(2,2) +  D_rot(3,3));  %timescale for loss of orientation around swimming axis, ignoring fact that this is not the principle axis and there are cross terms
    else
        inds = setdiff(1:3,varargin{1});
        tau_a(i) = 1/ (D_rot(inds(1),inds(1)) +  D_rot(inds(2),inds(2)));
    end
    
end