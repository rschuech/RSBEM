
AR2 = Metadata(1).geom.nturns/2;  r = Metadata(1).geom.radius_curv;  y = 2*r*(sin(AR2*pi))^2, x = 2*r*sin(AR2*pi)*cos(AR2*pi)
%from wikipedia chord length formula


center = NaN(3,length(timestepping_solution.x));
center0 = [x; y; 0];

%tic

parfor i = 1:length(timestepping_solution.x)
    
    rotmat =   rotation_matrix('z',timestepping_solution.y(i,6)) * rotation_matrix('y',timestepping_solution.y(i,5)) * rotation_matrix('x',timestepping_solution.y(i,4));
    
    center(:,i) = rotmat * center0  +  timestepping_solution.y(i,1:3)';
    
end


%toc

