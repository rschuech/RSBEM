function [pt] = tail_pt_wrapper(t,time,geom)

 [pt, vel, der] = tail_parameterized(t, time, geom);
 
[pt, vel, der] = change_ref_frame(pt, vel, der, geom);

