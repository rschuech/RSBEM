function [pt, vel, der] = change_ref_frame(pt, vel, der, geom)

for i = 1:size(pt, 1)
    
    pt(i,:) = [-1 0 0; 0 -1 0; 0 0 1] * pt(i,:)';% rotate 180 deg around z-axis as in script
    pt(i,:) = pt(i,:) + geom.translation';  %translate to final location as in script
    
    
    vel(i,:) = [-1 0 0; 0 -1 0; 0 0 1] * vel(i,:)';% rotate 180 deg around z-axis as in script
    
    
    der(i,:) = [-1 0 0; 0 -1 0; 0 0 1] * der(i,:)';% rotate 180 deg around z-axis as in script
    
end