figure

rng(23); nxy = 31;
xy = 2*(rand(2,nxy)-.5); vals = sum(xy.^2);
noisyvals = vals + (rand(size(vals))-.5)/5;
st = tpaps(xy,noisyvals); fnplt(st), hold on
avals = fnval(st,xy);
plot3(xy(1,:),xy(2,:),vals,'wo','markerfacecolor','k')
quiver3(xy(1,:),xy(2,:),avals,zeros(1,nxy),zeros(1,nxy), ...
    noisyvals-avals,'r'), hold off

%%
clear sweep

metric = 'eff';
vars = {'amp' , 'lambda' , 'nlambda'};
for v = 1:length(vars)
    var = vars{v};
    Points = F.(metric).F.(var).Points;
    Values = F.(metric).F.(var).Values;
    
    st = tpaps(Points',Values');
    
%     Points_interp = Points;
    Points_interp = standardize( [   1.4525 0;   3.875 0.7125;    3.875 0.725;   3.875 0.6875;   4.125 0.7375;  4.125 0.7125;  4.125 0.7;  4.125 0.6875] , limits) ;
    
    sweep.(var) = fnval(st,Points_interp')';
    
    figure(384);
    fnplt(st);  shading interp
%     pause
    
end
sweep.bodies = unstandardize( Points_interp , limits);

Bodies = sweep.bodies;  Tails = [sweep.amp sweep.lambda sweep.nlambda];
%%
racks = [1  2  3  4];
reps = [1 1 1 1];
% reps = [1 2 1 2 1 2 1 2];

for n = 1:length(racks)
    clear sweep
    % n = 4;
    bodies = Bodies(n:length(racks):end,:);
    guesses = Tails(n:length(racks):end,:);
%     inds = randperm(size(bodies,1));
%     bodies = bodies(inds,:);
%     guesses = guesses(inds,:);
sweep.bodies = bodies;  
sweep.amp = guesses(:,1); sweep.lambda = guesses(:,2);  sweep.nlambda = guesses(:,3);

    save(['C:\Hull\efficiency_smooth_refined_CFD0',num2str(racks(n)),'.mat'],'sweep');

    
end
