


function [sqdist] = dist(xi_eta,x0,element_nodes,shape_parameters)

x = T6interp(element_nodes,xi_eta(1),xi_eta(2),shape_parameters);

sqdist = sum( (x - x0).^2 );