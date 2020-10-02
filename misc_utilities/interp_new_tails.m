clear interp_input

interp_input.AR1_AR2 = [1 0;  10 0; ...
    1.1 0.35; 5 0.35; 10 0.35; ...
    2 0.635; ...
    5.75 0.65; 10 0.65; ...
    2.5 0.795; ...
    3 0.905; ...
    3.5 0.945; ...
    5.5 0.985; ...
    10 0.95; ...
    10 0.5;];

vars = {'amp','lambda','nlambda'};
AR1_AR2_Results = [ [Results.AR1]' [Results.AR2]' ];

for i = 1:size(interp_input.AR1_AR2,1)
    
    [~,ind] = ismember(interp_input.AR1_AR2(i,:),AR1_AR2_Results,'rows');
    
    
    for v = 1:length(vars)
        var = vars{v};
        interp_input.(var)(i) = Results(ind).(var);
    end
    
end


[AR1_AR2,inds] = unique([ [Results.AR1]' [Results.AR2]' ],'rows');
AR1_range = max(AR1_AR2(:,1)) - min(AR1_AR2(:,1));
AR2_range = max(AR1_AR2(:,2)) - min(AR1_AR2(:,2));
limits = [min(AR1_AR2(:,1)) max(AR1_AR2(:,1)); min(AR1_AR2(:,2)) max(AR1_AR2(:,2))];
[standardized] = standardize(interp_input.AR1_AR2, limits);

[standardized_XY] = standardize([X(:) Y(:)], limits);
standardized_X = NaN(size(X));  standardized_X(:) = standardized_XY(:,1);
standardized_Y = NaN(size(Y));  standardized_Y(:) = standardized_XY(:,2);

for v = 1:length(vars)
    var = vars{v};
    %     method = 'linear';
    method = 'natural';
    
    
    
    F = scatteredInterpolant(standardized(:,1), standardized(:,2), interp_input.(var)',method,'linear');
    interpolants.(var) = F;
    
    
%     Z = F(standardized_X,standardized_Y)  / F(0,0);
     Z = F(standardized_X,standardized_Y)  ;
    
    figure(500 + v)
    clf
    
    pcolor(X,Y,Z)
    %shading interp
    shading flat
    
    set(gca,'fontsize',fontsizes.axes);
    
    
    xlabel('AR_1 (elongation)','fontsize',fontsizes.labels)
    ylabel('AR_2 (curvature)','fontsize',fontsizes.labels)
    cbar = colorbar;
    set(gcf,'Position',[ 680          82        1051         896]);
    ylim([0 1]);
    set(gca,'XTick',[1:10])
    hold on
    d = plot(X_Y_unq(subset,1),X_Y_unq(subset,2), '.','markersize',3,'markerfacecolor','none','markeredgecolor','k');
    d2 = plot(interp_input.AR1_AR2(:,1),interp_input.AR1_AR2(:,2), 'o','markersize',8,'markerfacecolor','k','markeredgecolor','k');
    hold off
    
end


%%


AR1_AR2_new = AR1_AR2_Results;
AR1_AR2_new(AR1_AR2_new(:,1)==6.25 & AR1_AR2_new(:,2)==0.675 , :) = [];
AR1_AR2_new(AR1_AR2_new(:,1)==6.25 & AR1_AR2_new(:,2)==0.625 , :) = [];
AR1_AR2_new(end+1,:) = [5.5 0.625];
AR1_AR2_new(end+1,:) = [5.5 0.625];

[standardized_new] = standardize(AR1_AR2_new, limits);

tails_new.amp = interpolants.amp(standardized_new(:,1),standardized_new(:,2));
tails_new.lambda = interpolants.lambda(standardized_new(:,1),standardized_new(:,2));
tails_new.nlambda = repmat( interp_input.nlambda(1) , size(AR1_AR2_new,1),1);  % first entry of interp_input should be a sphere, which we normalize to and nlambda seems nearly constant around this value

bodies0 = AR1_AR2_new;
guesses0 = [tails_new.amp tails_new.lambda tails_new.nlambda];
%%
racks = [1  2  3 ];
reps = [1 1 1];
% reps = [1 2 1 2 1 2 1 2];
allbods = [];
for n = 1:length(racks)
    % n = 4;
    bodies = bodies0(n:length(racks):end,:);
    guesses = guesses0(n:length(racks):end,:);
    inds = randperm(size(bodies,1));
    bodies = bodies(inds,:);
    guesses = guesses(inds,:);
    save(['E:\Hull\sweep_CFD0',num2str(racks(n)),'_',num2str(reps(n)),'.mat'],'bodies','guesses');
    allbods = [allbods; bodies];
    
end