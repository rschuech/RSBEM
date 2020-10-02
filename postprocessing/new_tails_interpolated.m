% manually define F_amp, F_lambda, F_nlambda by running
% swimming_efficiency2.m shat a few times


[ia, ib] = ismember(X_Y_unq(:,1), [1.125 1.25 1.375 1.5 1.625 1.75]);  %current bodies to exclude

new_bodies = X_Y_unq;
new_bodies(ia,:) = [];
new_bodies(new_bodies(:,1) > 10 , :) = [];  %get rid of AR1 > 10 too

AR1 = [ repmat( 1.1:0.1:2 , length([0.05:0.05:1]), 1) ];  %refinement in low AR1 region
AR2 = repmat([0.05:0.05:1]', 1,length(1.1:0.1:2));
bods = [ [AR1(:) AR2(:)] ;  [(1.1:0.1:2)' repmat(0.375,length(1.1:0.1:2),1)] ; [(6.5:0.5:10)' repmat(0.675,8, 1)]; [(6.5:0.5:10)' repmat(0.625,8, 1)]; ];
[~, obj_signed] = curved_rod_fat_limit(bods(:,1), bods(:,2), 1);
bods = bods(obj_signed > 0, :);

new_bodies = unique([new_bodies;   bods;],'rows');

%% OR, hardcode new_bodies here
AR2 = 0:0.05:0.35;
AR1 = repmat(1.25,size(AR2));
new_bodies = [ AR1' AR2'];
AR2 = 0:0.05:0.55;
AR1 = repmat(1.75,size(AR2));
new_bodies = [new_bodies; [AR1' AR2'] ];
AR1 = 1:0.25/4:2;  AR1 = setdiff(AR1,[1:0.25:2]);
AR2 = zeros(size(AR1));
new_bodies = [new_bodies; [AR1' AR2'] ];
%% Don't bother with this, use cell called
% interpolate tail guesses for secondary sweep from primary sweep results
% in swimming_efficiency.m

 % run swimming_efficiency for efficiency here to generate refined X_Y_unq
[AR1, AR2] = ndgrid( 1:0.5:10, 0:0.05:1);
AR1_AR2 = roundn([AR1(:) AR2(:)],-10);

[AR1, AR2] = ndgrid( 1:1:10, 0:0.1:1);
AR1_AR2_coarse = roundn([AR1(:) AR2(:)],-10);
% V = 1;
% [obj, obj_signed] = curved_rod_fat_limit(AR1_AR2(:,1), AR1_AR2(:,2), V);
% AR1_AR2 = AR1_AR2(obj_signed >= 0);  % remove self intersecting fat donuts
new_bodies = setdiff(  intersect(AR1_AR2,X_Y_unq,'rows')  ,  AR1_AR2_coarse ,'rows' );

%%
limits = [min(X(:)) max(X(:)); min(Y(:)) max(Y(:))];
if max(limits(:)) > 10  %make sure limits are right for new AR1 cutoff
    stopafro
end

found = 0;  notfound = 0;
clear bodies guesses
for i = 1:size(new_bodies,1)  % data points we've loaded and plotted
    %     i
    body = new_bodies(i,:);
    %     body_std = [standardized_X_subset(i) standardized_Y_subset(i)];
    body_std = standardize(body,limits);
    
    
    Fs = {F_amp F_lambda F_nlambda};
    clear new_tail;
    for j = 1:length(Fs)
        
        F = Fs{j};
        ind = find(ismember(F.Points,body_std,'rows'));
        if numel(ind) < 1
            notfound = notfound + 1;
             new_tail(j) = F(body_std);
            
        elseif numel(ind) > 1
            stopafra
        else
            
            % otherwise, ind == 1 and we found exactly one match for this point in
            % F
            %delete current point from interpolant so it's effect is ignored, in
            %case it's an outlier
            F.Points(ind,:) = [];
            F.Values(ind) = [];
            found = found + 1;
            new_tail(j) = F(body_std);
        end
        
    end
    
    %     old_tail = [best_amps(i) best_lambdas(i) best_nlambdas(i)];
    
    
    guesses(i,:) = new_tail;
    bodies(i,:) = body;
    
    
end


clear new_sweep
new_sweep.AR1 = bodies(:,1)';
new_sweep.AR2 = bodies(:,2)';
new_sweep.amp = guesses(:,1)';
new_sweep.lambda = guesses(:,2)';
new_sweep.nlambda = guesses(:,3)';
