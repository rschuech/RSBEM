
% load('C:\Users\rudi\Desktop\RD\Results\Results_Tumbling_Ease_SDE.mat');
folder = 'C:\Users\rudi\Desktop\RD\Tumbling_dumps_all\';
% folder = 'Z:\Tumbling_dumps_all\';


files = dir([folder,'*.mat']);  files = {files.name};
Results = [];

for di = 255:338 % 1:84  85:169  170:254
    di
    file = files{di};
    load([folder, file],'D','avg_swimming_axes','input');
    
    Dm = D;
    body = input.body.AR;
    
    %  D = [0.5 0.01 0.05]';
    % D = [0.5 0.25 1.5]';
    D = diag(Dm.rotation.diffusivity);
    
    % principal_axes = [ [1 0 0]' [0 1 0]' [0 0 1]' ];
    principal_axes = Dm.rotation.axes;
    % principal_axes(:,2) = - principal_axes(:,2);
    
    % swimming_axis0 = [1 0 0]';
    % swimming_axis0 = principal_axes(:,3);
    swimming_axis0 = avg_swimming_axes.avg_swimming_axis(avg_swimming_axes.AR1_AR2(:,1) == body(1) & avg_swimming_axes.AR1_AR2(:,2) == body(2) , :)';
    swimming_axis0 = swimming_axis0 / sqrt(sum(swimming_axis0.^2));
    
    dt = fnval(pp_dt,body(1));
    
    
    
    % Nt = 1E4;
    
    Tmax = fnval(pp_T,body(1));
    T = 0:dt:Tmax;
    Nt = length(T);
    
    Np = 2E7;
    % dx = randn(Nt,Np, length(D)) .* repmat( reshape( sqrt(2*D*dt) , [1 1 length(D)] ) , Nt,Np );  % can either be change in distance or theta
    
    
    X0 = [0 0 0]';
    
    [    cos_theta ] = rotational_diffusion_mex(Nt,Np,dt,principal_axes, swimming_axis0,D);
    
    
    %
    %%
    mean_cos_theta = mean(cos_theta,2);
    % [fitobject,gof] = fit(T(1:end)',mean_cos_theta(1:end),'exp1','Lower',[1 -Inf],'Upper',[1 Inf]);
    % [fitobject,gof] = fit(T(1:end)',mean_cos_theta(1:end),'exp1','Lower',[1 -0.6135],'Upper',[1 -0.6135]);
    figure(66);
    plot(T,mean_cos_theta,'o-'); grid on;
    hold on;
    % plot(T,fitobject.a*exp(fitobject.b*T),'k--','linewidth',2);
    pp = csaps(T,mean_cos_theta);
    plot(T,fnval(pp,T),'k--','linewidth',2);
    hold off
    
    
    % tau = -1/fitobject.b
    
    tau = fzero( @(t) fnval(pp,t) - 1/exp(1) , 5);
    title(num2str(tau));
    drawnow
    
    assert(tau < 0.95*T(end),'tau near T(end)');
    assert(tau > 0,'tau < 0)');
    
    
    Results(di).AR1 = body(1);  Results(di).AR2 = body(2);
    Results(di).amp = NaN;  Results(di).lambda = NaN;  Results(di).nlambda = NaN;
    Results(di).rotational_diffusivity = Dm.rotation.diffusivity;
    Results(di).principle_axis_1 = Dm.rotation.axes(:,1);
    Results(di).avg_swimming_axis = swimming_axis0;
    Results(di).tau_a = tau;
    
    if mod(di,5) == 0
        save('C:\Users\rudi\Desktop\RD\Results\Results_Tumbling_Ease_SDE.mat','Results');
    end
    
end
 save('C:\Users\rudi\Desktop\RD\Results\Results_Tumbling_Ease_SDE.mat','Results');

return

%%

L = [1 2 4 6 8 10];
dt =[ 0.1 0.1 0.15 0.25 0.3 0.4];
T = [1.5 1.5 3 5 8 10];
L_refined = linspace(1,10,100);

pp_dt = pchip(L,dt);
pp_T = pchip(L,T);

figure(467)
plot(L,dt,'bo','markerfacecolor','b'); hold on;
plot(L,T,'ro','markerfacecolor','r');
plot(L_refined,fnval(pp_dt,L_refined),'b-');
plot(L_refined,fnval(pp_T,L_refined),'r-');
grid on;  hold off;
