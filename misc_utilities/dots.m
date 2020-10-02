function [dots] = dots(in)

% takes dot products of all three vectors of in matrix with each other

dots(:,1) = dot(in(:,1),in(:,2));
dots(:,2) = dot(in(:,1),in(:,3));
dots(:,3) = dot(in(:,2),in(:,3));

