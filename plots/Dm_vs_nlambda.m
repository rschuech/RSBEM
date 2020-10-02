


tails = [ [Results.amp]'  [Results.lambda]'  [Results.nlambda]'];
Dms = [Results.Dm];

lambda_3 = find(tails(:,1) == 0.4 & tails(:,2) == 3);
lambda_4 = find(tails(:,1) == 0.4 & tails(:,2) == 4);

figure(45);  clf
plot(tails(lambda_3,3),Dms(lambda_3),'ro-','markerfacecolor','r');
hold on
plot(tails(lambda_4,3),Dms(lambda_4),'bo-','markerfacecolor','b');
grid on
legend('amp = 0.4  lambda = 3','amp = 0.4  lambda = 4','location','best');
xlabel('nlambda')
ylabel('Dm')

%%
results = [[Results.AR1]' [Results.AR2]' [Results.amp]'  [Results.lambda]'  [Results.nlambda]'];
Dms = [Results.Dm];

AR5 = find(results(:,1) == 5);


figure(46);  clf
plot(results(AR5,end),Dms(AR5),'ro-','markerfacecolor','r');

grid on
%legend('amp = 0.4  lambda = 3','amp = 0.4  lambda = 4','location','best');
xlabel('nlambda')
ylabel('Dm')

%%
results = [[Results.AR1]' [Results.AR2]' [Results.amp]'  [Results.lambda]'  [Results.nlambda]'];
Dms = [Results.Dm];

AR5 = find(results(:,1) == 5 & results(:,2) == 0.65);


figure(47);  clf
plot(results(AR5,end),Dms(AR5),'ro-','markerfacecolor','r');

grid on
%legend('amp = 0.4  lambda = 3','amp = 0.4  lambda = 4','location','best');
xlabel('nlambda')
ylabel('Dm')


%%

results = [[Results.AR1]' [Results.AR2]' [Results.amp]'  [Results.lambda]'  [Results.nlambda]'];
Dms = [Results.Dm];

amps = find(results(:,1) == 1 & results(:,2) == 0 & results(:,4) == 3.5 & results(:,5) == 1.3);
x = results(amps,3);  y = Dms(amps);
[x,inds] = sort(x);  y = y(inds);

figure(51);  clf
plot(x,y,'ro-','markerfacecolor','r');

grid on
%legend('amp = 0.4  lambda = 3','amp = 0.4  lambda = 4','location','best');
xlabel('amp')
ylabel('Dm')


%%

results = [[Results.AR1]' [Results.AR2]' [Results.amp]'  [Results.lambda]'  [Results.nlambda]'];
Dms = [Results.Dm];

lambdas = find(results(:,1) == 1 & results(:,2) == 0 & results(:,3) == 0.5 & results(:,5) == 1.3);
x = results(lambdas,4);  y = Dms(lambdas);
[x,inds] = sort(x);  y = y(inds);

figure(50);  clf
plot(x,y,'ro-','markerfacecolor','r');

grid on
%legend('amp = 0.4  lambda = 3','amp = 0.4  lambda = 4','location','best');
xlabel('lambda')
ylabel('Dm')