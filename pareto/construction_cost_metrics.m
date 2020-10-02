%% add Constrution Ease to pareto_data struct
limits = [1 10; 0 1];  npts = 300;
method = 'natural';
% method = 'linear';
[X,Y] = meshgrid(linspace(limits(1,1),limits(1,2),npts),linspace(limits(2,1),limits(2,2),npts));
[standardized_XY] = standardize([X(:) Y(:)], limits);
standardized_X = NaN(size(X));  standardized_X(:) = standardized_XY(:,1);
standardized_Y = NaN(size(Y));  standardized_Y(:) = standardized_XY(:,2);


[curvs] = curved_rod_curvature(X, Y, 1);
[area, radius_curv1, radius_curv2] = curved_rod_area(X,Y,1); % curv1 is circumferential (depends on SF1), curv2 is axial (depends on SF2 but also SF1 due to skinnier donuts having less curvature than fat donuts)
% radius_curv1(isinf(radius_curv1)) = max(radius_curv1(~isinf(radius_curv1))) * 1.01;  
curvs.curv1 = 1./radius_curv1;  % replace Inf by arbitrary large number, and invert to get a curvature like rest of metrics
% radius_curv2(isinf(radius_curv2)) = max(radius_curv2(~isinf(radius_curv2))) * 1.01;  
curvs.curv2 = 1./radius_curv2;  % replace Inf by arbitrary large number, and invert to get a curvature like rest of metrics

  V = 1;  %volume, microns^3
 sph_rad = (V*3/4/pi)^(1/3);


 
clear Construction_Eases
%% unweighted and weighted main curvature cost options
fields = fieldnames(curvs);
for f = 1:length(fields)
    % unweighted
%     Z = 1./curvs.(feields{f});  
%     Z = - curvs.(fields{f}) + max(curvs.(fields{f})(:));
    if any(0 == curvs.(fields{f})(:))  % this curvature metric is sometimes zero, need to prevent blowup of 1/curv
         C = 1/sph_rad;  % fudge constant to add to prevent 1/0 blowups
    else
        C = 0;
    end
%    1 * exp( - 1 * curvs.(fields{f}) )  
   Z = 1 ./ ( curvs.(fields{f}) + C ) ;
%        Z = Z/max(Z(:));
    F = scatteredInterpolant(standardized_X(~isnan(Z)), standardized_Y(~isnan(Z)), Z(~isnan(Z)),method,'none');
    Construction_Eases.unweighted.(fields{f}).F.metric = F;
    
    % weighted, cost is multiplied by surface area, then the whole thing
    % inverted for Ease
%     Z = 1./ ( curvs.(fields{f}) .* area );  
%     Z = - curvs.(fields{f}) .* area;  Z = Z + max(-Z(:));
    Z = sqrt(   1 ./ ( curvs.(fields{f}) + C ) ./ area   );  %
    
%     Z = round(Z,3,'significant');
%     Z = Z/max(Z(:));
    F = scatteredInterpolant(standardized_X(~isnan(Z)), standardized_Y(~isnan(Z)), Z(~isnan(Z)),method,'none');
    Construction_Eases.weighted.(  fields{f}  ).F.metric = F;
end

% construction ease based on inverse surface area alone
Z = 1./area;    % Z = Z/max(Z(:));
%  Z = round(Z,3,'significant');
F = scatteredInterpolant(standardized_X(~isnan(Z)), standardized_Y(~isnan(Z)), Z(~isnan(Z)),method,'none');
Construction_Eases.area_alone.F.metric = F;


%  Z = 1 ./ ( curvs.curv1 .* curvs.curv2  +1) ./ area;
  Z = 1 ./ ( curvs.curv1 .* curvs.curv2  + C) ./ area ;
 
%   Z = round(Z,3,'significant');
  
 F = scatteredInterpolant(standardized_X(~isnan(Z)), standardized_Y(~isnan(Z)), Z(~isnan(Z)),method,'none');
    Construction_Eases.weighted.test.F.metric = F;
%%

save('C:\Hull\Results\fixed\Construction_Ease_fixed_scaling.mat','Construction_Eases');



