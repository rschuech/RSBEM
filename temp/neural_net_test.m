%%

x = [eff.curved_rod.AR1' eff.curved_rod.AR2' eff.curved_rod.amp' eff.curved_rod.lambda' eff.curved_rod.nlambda'];
y = eff.curved_rod.power_effs;

pm = polyfitn(x, y, 3);

abserrors = abs( (polyvaln(pm,x) - y') ./ y');
[min(abserrors) mean(abserrors) max(abserrors)]*100

%%

range.AR1 = Inf;  %use data + and - this amount on either side of current AR1
range.AR2 = Inf;

AR1s = [2.7 2.9   3     4.75   4.75   4.75    5.25   5.25   5.25];
AR2s = [0   0    0.85   0.625  0.65   0.675   0.625  0.65  0.675];

% AR1s = [];  AR2s = [];

[AR1_AR2,inds] = unique([ [Results.AR1]' [Results.AR2]' ; AR1s' AR2s'],'rows');

pms_AR = AR1_AR2;

AR1_range = max(AR1_AR2(:,1)) - min(AR1_AR2(:,1));
AR2_range = max(AR1_AR2(:,2)) - min(AR1_AR2(:,2));
 
    inc.AR1 = 0.25;  inc.AR2 = 0.025;
clear pms
errors = NaN(size(AR1_AR2,1),3);
lastwarn('');
parfor i = 1:size(AR1_AR2,1) %each unique body
    if rem(i,20) == 0
        i/size(AR1_AR2,1)
    end
    body = AR1_AR2(i,:);

    inds = [Results.AR1] >= body(1) - range.AR1 & [Results.AR1] <= body(1) + range.AR1        &         [Results.AR2] >= body(2) - range.AR2 & [Results.AR2] <= body(2) + range.AR2;
    
    %inds = best_inds;  %single best tails for all bodies
    
    
%     x = [[Results(inds).AR1]' [Results(inds).AR2]' [Results(inds).amp]' [Results(inds).lambda]' [Results(inds).nlambda]'];
%     y = [Results(inds).Power_eff]';
    
    x = [X_Y_unq  best_amps best_lambdas best_nlambdas];
    y = best_effs;
    
    standardized_dists = sqrt(  1/AR1_range^2*( x(:,1) - body(1) ).^2  +  1/AR2_range^2*( x(:,2) - body(2) ).^2  );
    ep = 20;  %bigger this is, sharper rbf curve drops off away from current body point
    % ep = 0.2;   % was 100
    % ep = 10;
     
     % ep = 5;
     
    rbf = @(r) 1./(1+(ep*r).^2); %see wikipedia on radial basis functions
    weights = rbf(standardized_dists);
   
    pms(i) = polyfitn(x, y, 2, weights);
    
  
%     c = 0;
%     while ~isempty(lastwarn)
%         lastwarn('');  %innocent until proven quilty
%         c = c + 1
%         inds = [Results.AR1] >= body(1) - range.AR1 - c*inc.AR1 & [Results.AR1] <= body(1) + range.AR1 + c*inc.AR1       &         [Results.AR2] >= body(2) - range.AR2 - inc.AR2 & [Results.AR2] <= body(2) + range.AR2 + c*inc.AR2;
%         
%         x = [[Results(inds).AR1]' [Results(inds).AR2]' [Results(inds).amp]' [Results(inds).lambda]' [Results(inds).nlambda]'];
%         y = [Results(inds).Power_eff]';
%         
%     stopooo
%         pms(i) = polyfitn(x, y, 2);
%     end
    
    abserrors = abs( (polyvaln(pms(i),x) - y) ./ y);
  errors(i,:) =   [min(abserrors) mean(abserrors) max(abserrors)]*100;
    if length(y) <= 3
        %     pause
    end
    
end





%%
N = 20;


spread = 1;

x = [eff.curved_rod.AR1; eff.curved_rod.AR2; eff.curved_rod.amp; eff.curved_rod.lambda; eff.curved_rod.nlambda;];
y = eff.curved_rod.power_effs;

tic;  net_0p025 = newrbe(x,y,spread); toc

abserrors = abs( ( net(x) - y) ./ y);
[min(abserrors) mean(abserrors) max(abserrors)]*100




tic;  netcc = newrbe(x,y,4); toc

tic;  netccc = newrbe(x,y,8); toc

tic;  net = newrbe(x,y,1); toc



tic;  netrr = newrbe(x,y,0.25); toc

%st = tpaps(x(:,1:N),y(1:N))

sc = 1 * 1;    % spread constant

tic;  net = newrbe(x,y,1); toc


sc = 1 * 2;    % spread constant

tic;  netc = newrbe(x,y,sc); toc

tic;  netcc = newrbe(x,y,4); toc


%der = defaultderiv('dperf_dwb',net,x(:,2),y(2));

net2 = newrbe(x,y,0.5);
% 3.2664e-08     0.067789       5.4678

net3 = newrb(x,y,1E-6, 1);

net4 = newrbe(x,y,0.25);
%   1.0742e-08     0.024203       2.8852

net5 = newrbe(x,y,0.21);

abserrors = abs( ( netr(x) - y) ./ y);
[min(abserrors) mean(abserrors) max(abserrors)]*100
%%
err_cutoff = 9/100;

[eff.curved_rod.AR1(abserrors >= err_cutoff)' eff.curved_rod.AR2(abserrors >= err_cutoff)' eff.curved_rod.amp(abserrors >= err_cutoff)' eff.curved_rod.lambda(abserrors >= err_cutoff)' eff.curved_rod.nlambda(abserrors >= err_cutoff)']

%%


obj = @(x) net4(x);
% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
solution2 = fmincon(obj, [5 0.6 0.4 2 1.5]',[],[],[],[],[1 0 0 0 0]',[12 1 10 10 10]')


%%
AR1 = 1.5;  AR2 = 0.1;

obj = @(x) -net([AR1; AR2; x]);

lb = [ min(eff.curved_rod.amp) min(eff.curved_rod.lambda) min(eff.curved_rod.nlambda)]';
ub = [max(eff.curved_rod.amp) max(eff.curved_rod.lambda) max(eff.curved_rod.nlambda)]';

lb = [ min(eff.curved_rod.amp) min(eff.curved_rod.lambda) 1]';
ub = [0.4 max(eff.curved_rod.lambda) max(eff.curved_rod.nlambda)]';

tail = fmincon(obj, [ 0.4 2 1.5]',[],[],[],[],lb,ub)


%%
% lambda = 2.177;
lambda = 2.185;
nlambda = 1.863;

amp = linspace(0.02,1.5,500);

xx = [ repmat(AR1,1,length(amp)); repmat(AR2,1,length(amp)); amp; repmat(lambda,1,length(amp)); repmat(nlambda,1,length(amp));];

ef = netr(xx) / max(eff.curved_rod.power_effs);

hold on
pl = plot(amp,ef,'b--');
hold off

%%

amp = plotdata.amp;

xx = [ repmat(AR1,1,length(amp)); repmat(AR2,1,length(amp)); amp; repmat(lambda,1,length(amp)); repmat(nlambda,1,length(amp));];

ef = net(xx) / max(eff.curved_rod.power_effs);


abserrors = abs( ( ef - plotdata.amp_deps) ./ plotdata.amp_deps) * 100











%% guess first tail for brand new bodies




range.AR1 = Inf;  %use data + and - this amount on either side of current AR1
range.AR2 = Inf;

[AR1_AR2,inds] = unique([ [Results.AR1]' [Results.AR2]' ],'rows');

AR1_range = max(AR1_AR2(:,1)) - min(AR1_AR2(:,1));
AR2_range = max(AR1_AR2(:,2)) - min(AR1_AR2(:,2));


AR1s = [1.025    1.01     1.2     1.975     2.5    3       3.75     9.5     10.5      7.25      11.5       2.7      2.3      10.5     11.5     11.5];
AR2s = [0.055    0.2      0.375   0.625     0.25   0.825   0.95     0.85     0.9      0.94      0.9        0.85     0.65     0.05     0.05     0.75];


AR1_AR2 = [AR1s' AR2s'];

clear pms_new
errors = NaN(size(AR1_AR2,1),3);
lastwarn('');
for i = 1:size(AR1_AR2,1) %each unique body
    if rem(i,1) == 0
        i/size(AR1_AR2,1)
    end
    
    body = AR1_AR2(i,:);

    inds = [Results.AR1] >= body(1) - range.AR1 & [Results.AR1] <= body(1) + range.AR1        &         [Results.AR2] >= body(2) - range.AR2 & [Results.AR2] <= body(2) + range.AR2;
    
    %inds = best_inds;  %single best tails for all bodies
    
    
    x = [[Results(inds).AR1]' [Results(inds).AR2]' [Results(inds).amp]' [Results(inds).lambda]' [Results(inds).nlambda]'];
    y = [Results(inds).Power_eff]';
    
    standardized_dists = sqrt(  1/AR1_range^2*( x(:,1) - body(1) ).^2  +  1/AR2_range^2*( x(:,2) - body(2) ).^2  );
    ep = 20;  %bigger this is, sharper rbf curve drops off away from current body point
     ep = 0.2;   % was 100

     ep = 20;
     
    rbf = @(r) 1./(1+(ep*r).^2); %see wikipedia on radial basis functions
    weights = rbf(standardized_dists);
   
    pms_new(i) = polyfitn(x, y, 2, weights);
    
    
    inc.AR1 = 0.25;  inc.AR2 = 0.025;
    c = 0;
    while ~isempty(lastwarn)
        lastwarn('');  %innocent until proven quilty
        c = c + 1
        inds = [Results.AR1] >= body(1) - range.AR1 - c*inc.AR1 & [Results.AR1] <= body(1) + range.AR1 + c*inc.AR1       &         [Results.AR2] >= body(2) - range.AR2 - inc.AR2 & [Results.AR2] <= body(2) + range.AR2 + c*inc.AR2;
        
        x = [[Results(inds).AR1]' [Results(inds).AR2]' [Results(inds).amp]' [Results(inds).lambda]' [Results(inds).nlambda]'];
        y = [Results(inds).Power_eff]';
        
        pms_new(i) = polyfitn(x, y, 2);
    end
    
    abserrors = abs( (polyvaln(pms_new(i),x) - y) ./ y);
  errors(i,:) =   [min(abserrors) mean(abserrors) max(abserrors)]*100;
    if length(y) <= 3
        %     pause
    end
    
end


%% now sample polynomials for best tail guesses
    lb = [ min([Results.amp]) min([Results.lambda]) min([Results.nlambda])]';
    ub = [max([Results.amp]) max([Results.lambda]) max([Results.nlambda])]';
    
    clear next_iter
    
for i = 1:size(AR1_AR2,1) %each unique body

    
    body = AR1_AR2(i,:);
    
    
 pm = pms_new(i);
     obj2 = @(x)  - polyvaln(pm,[body(1) body(2) x']);
     
     pm_guess = fmincon(obj2, [1 1 1]',[],[],[],[],lb,ub);
    
    next_iter(i).param = 'all';
%     next_iter(body_i).nlambda = net_guess(3);
%     next_iter(body_i).amp = net_guess(1);
%     next_iter(body_i).lambda = net_guess(2);
    
        next_iter(i).nlambda = pm_guess(3);
    next_iter(i).amp = pm_guess(1);
    next_iter(i).lambda = pm_guess(2);
    
end

%%  add on these brand new body/tails onto end of new_sweep

new_sweep.AR1 = [new_sweep.AR1 AR1s];
new_sweep.AR2 = [new_sweep.AR2 AR2s];
new_sweep.amp = [new_sweep.amp  [next_iter.amp] ];
new_sweep.lambda = [new_sweep.lambda  [next_iter.lambda] ];
new_sweep.nlambda = [new_sweep.nlambda  [next_iter.nlambda] ];

