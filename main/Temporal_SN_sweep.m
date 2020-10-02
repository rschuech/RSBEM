rack = getenv('computername');



num = '1';

sweep_file = ['sweep_',rack,'_',num2str(num),'.mat'];



AR1 = 1;  AR2 = 0;


amps = linspace(0.05,1,10);
lambdas = linspace(0.5,5,10);
nlambdas = linspace(0.1,3,10);

[Amp,Lambda,Nlambda] = ndgrid(amps,lambdas,nlambdas);
Arclength = NaN(size(Amp));
for i = 1:numel(Amp)
    Arclength(i) = bacterial_tail_arclength([Amp(i) Lambda(i) Nlambda(i)]);
end

filter = (Arclength(:) > 5 & Arclength(:) < 10);

tails = [Amp(filter) Lambda(filter) Nlambda(filter)];
arclengths = Arclength(filter);
%%

lock_file_name = ['results_lock_',num];
results_file_name =   ['Results_Temporal_SN_sweep_',num,'.mat'];

% results_file_name = 'Results_.mat';

sweep_tempfile = ['tempsweep_',num2str(num),'.mat'];

results_file = ['C:\Users\rudi\Desktop\RD\Results\',results_file_name];

Temporal_SN = NaN(size(tails,1),1);

%%
for sw = 1:size(tails,1)     %25:71
    sw
    % AR1 = bodies_sweep(sw,1);  AR2 = bodies_sweep(sw,2);
    
    % disp(['On opt sweep iter ',num2str(sw),' out of ',num2str(size(bodies_sweep,1)),'      body AR1 ',num2str(AR1),'  AR2 ',num2str(AR2)]);
    tail = tails(sw,:);
    

Temporal_SN(sw) =  optimization_wrapper_temporal_SN(AR1, AR2, tail(1), tail(2), tail(3), results_file, lock_file_name, results_file_name , sweep_tempfile);
    
    arclengths_temp = arclengths(~isnan(Temporal_SN));
    tails_temp = tails(~isnan(Temporal_SN),:);
    SN_temp = Temporal_SN(~isnan(Temporal_SN));
    [~,inds] = sort(SN_temp,'descend');
    
    [SN_temp(inds)  arclengths_temp(inds)  tails_temp(inds,:)]
    
    results = [SN_temp(inds)  arclengths_temp(inds)  tails_temp(inds,:)];
    
end
%%



amps = linspace(0.05,1,50);
lambdas = linspace(0.5,5,50);
nlambdas = linspace(0.1,3,50);

[Amp,Lambda,Nlambda] = meshgrid(amps,lambdas,nlambdas);

interpolant = scatteredInterpolant(results(:,3),results(:,4),results(:,5),results(:,1),'natural','none');
interped = interpolant(Amp,Lambda,Nlambda);

Arclength = NaN(size(Amp));
for i = 1:numel(Amp)
    Arclength(i) = bacterial_tail_arclength([Amp(i) Lambda(i) Nlambda(i)]);
end




%% trying various choices of max arclength gives best Temporal SN at an arclength very close to the max, suggesting best tail is always max arclength tail
filter = Arclength <= 100;
temp = interped(filter);
[~,ind] = max(temp);
temp = Arclength(filter);
temp(ind)

%%
figure(824)
clf
isosurface(Amp,Lambda,Nlambda,interped,3000);
grid on
xlabel('amp'); ylabel('lambda'); zlabel('nlambda')
