% run precompute_Pareto.m with usual resolution, then run this




Opt = Optimal;  Opt(isnan(Opt)) = 0; % need binary 1, 0 array

 B = bwboundaries(Opt);
 
main_region = B{1};  island_region = B{2};

%% main
NE_corner_ind = find(main_region(:,1) == size(Opt,1),1,'first');

upper_boundary = [ xp(main_region(1:NE_corner_ind,1))'    yp(main_region(1:NE_corner_ind,2))' ];
mid_ind = find(upper_boundary(:,1) > 5,1,'first');  % index near SF1 = 5 on upper boundary
early_ind = find(upper_boundary(:,1) > 1.5,1,'first');  % index near SF1 = 5 on upper boundary
upper_boundary = standardize(upper_boundary,limits);
t = 1:size(upper_boundary,1);  % rough attempt at parameterizing the boundary points
t_refined = linspace(1,t(end),5000);

n_knots = 5;
options = slmset('interiorknots','free','knots',n_knots,'xy',[1 0]);
slm(1) = slmengine(t,upper_boundary(:,1),options);
% options = slmset('interiorknots','free','knots',n_knots,'xy',[1 0.00],'decreasing',[t(mid_ind) t(end)],'increasing',[t(1) t(early_ind)]);
options = slmset('interiorknots','free','knots',n_knots,'xy',[1 0.00; ],'decreasing',[t(mid_ind) t(end)]);

slm(2) = slmengine(t,upper_boundary(:,2),options);
clear upper_boundary_splined
upper_boundary_splined(:,1) = slmeval(t_refined,slm(1));
upper_boundary_splined(:,2) = slmeval(t_refined,slm(2));
upper_boundary_splined = unstandardize(upper_boundary_splined,limits);

main_poly = polyshape([upper_boundary_splined; 10 0; 1 0;]);
% upper_boundary_splined = unstandardize(upper_boundary_splined,limits);
try, delete(smoothed_main), end;
smoothed_main = plot(upper_boundary_splined(:,1),upper_boundary_splined(:,2),'k-','LineWidth',1);
%%
n_shapes = 8;
% t_shapes = linspace(min(t),max(t),n_shapes);
% xy_shapes = [slmeval(t_shapes,slm(1))'  slmeval(t_shapes,slm(2))'];
% try, delete(xysh), end;
% xysh = plot(xy_shapes(:,1),xy_shapes(:,2),'bo','MarkerFaceColor','b','MarkerSize',6);

arclength = @(t1,t2) integral( @(t) sqrt( slmeval(t,slm(1),1).^2 + slmeval(t,slm(2),1).^2 ) , t1, t2);

tot_arclength = arclength(min(t),max(t));
arclength_shapes = linspace(0,tot_arclength,n_shapes);
t_guesses = linspace(min(t),max(t),n_shapes);
clear t_shapes;  t_shapes(1) = min(t);
for n = 2:n_shapes
    objfun = @(t_sol) arclength(t_shapes(1),t_sol) - arclength_shapes(n) ;
    t_shapes(n) = fzero(objfun , t_guesses(n));
end

xy_shapes = unstandardize( [slmeval(t_shapes,slm(1))'  slmeval(t_shapes,slm(2))'] , limits);
try, delete(xysh), end;
% xysh = plot(xy_shapes(:,1),xy_shapes(:,2),'bo','MarkerFaceColor','b','MarkerSize',6);
%% island

island_boundary = standardize( [ xp(island_region(:,1))'    yp(island_region(:,2))' ] , limits);
t = 1:size(island_boundary,1);  % rough attempt at parameterizing the boundary points
t_refined = linspace(1,t(end),1000);

options = slmset('interiorknots','free','knots',4,'endconditions','periodic');

slm(1) = slmengine(t,island_boundary(:,1),options);
slm(2) = slmengine(t,island_boundary(:,2),options);
clear island_boundary_splined
island_boundary_splined(:,1) = slmeval(t_refined,slm(1));
island_boundary_splined(:,2) = slmeval(t_refined,slm(2));
island_boundary_splined = unstandardize(island_boundary_splined , limits);

island_poly = polyshape(island_boundary_splined);

% island_boundary_splined = unstandardize(island_boundary_splined,limits);
try, delete(smoothed_island), end;
smoothed_island = plot(island_boundary_splined(:,1),island_boundary_splined(:,2),'k-','LineWidth',1);

temp = polyshape(island_boundary_splined(:,1),island_boundary_splined(:,2));
[centroid_x,centroid_y] = centroid(temp);  % copied, hardcoded into parameter_space_gridded