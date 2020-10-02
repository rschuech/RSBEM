
% load('C:\Users\rudi\Desktop\RD\Results\Results_Tumbling_Ease_SDE.mat');


tic

load('C:\Users\rudi\Desktop\RD\Tumbling_dumps_all\curved_rod_AR1_8_AR2_0_forced_dump.mat','D','avg_swimming_axes');
Dm = D;


%  D = [0.5 0.01 0.05]';
D = [0.5 0.25 1.5]';
% D = diag(Dm.rotation.diffusivity);

principal_axes = [ [1 0 0]' [0 1 0]' [0 0 1]' ];
% principal_axes = Dm.rotation.axes;
% principal_axes(:,2) = - principal_axes(:,2);

swimming_axis0 = [1 0 0]';
% swimming_axis0 = principal_axes(:,3);  
% swimming_axis0 = avg_swimming_axes.avg_swimming_axis(avg_swimming_axes.AR1_AR2(:,1) == 5 & avg_swimming_axes.AR1_AR2(:,2) == 0.6 , :)';
swimming_axis0 = swimming_axis0 / sqrt(sum(swimming_axis0.^2));

dt = 0.1;



% Nt = 1E4;

Tmax = 2;
T = 0:dt:Tmax;
Nt = length(T);

Np = 2E7;
% dx = randn(Nt,Np, length(D)) .* repmat( reshape( sqrt(2*D*dt) , [1 1 length(D)] ) , Nt,Np );  % can either be change in distance or theta


X0 = [0 0 0]';

[    cos_theta ] = rotational_diffusion_mex(Nt,Np,dt,principal_axes, swimming_axis0,D);


% do_progmon = true;
% 
% cos_theta = NaN(Nt,Np);
% if do_progmon
% ppm = ParforProgressStarter2('shat', Np);
% end
% 
% parfor j = 1:Np % particle
%     current_axes = principal_axes;
%     swimming_axis = swimming_axis0;
%     cos_theta_temp = NaN(Nt,1); cos_theta_temp(1) = 1;
%     for i = 2:Nt % time
%         %     i/size(dx,1)
%         dx = randn(3,1) .* sqrt(2*D*dt);
%         
%         order = randperm(length(D));
% % order = [1 2 3];
%           for k = 1:length(order)
% %         rotmat = rotate_arbitrary_vector( current_axes(:,order(1)), dx(order(1))) * rotate_arbitrary_vector( current_axes(:,order(2)), dx(order(2))) * rotate_arbitrary_vector( current_axes(:,order(3)), dx(order(3)));
%        
%          rotmat = rotate_arbitrary_vector( current_axes(:,order(k)), dx(order(k)));
%         
%         
%         current_axes = rotmat * current_axes;   current_axes = current_axes ./ sqrt(sum(current_axes.^2));
%         swimming_axis = rotmat * swimming_axis;   swimming_axis = swimming_axis ./ sqrt(sum(swimming_axis.^2));
%         
%           end
%           
%           
% %         figure(612)
% %         quiver3(0,0,0,swimming_axis(1),swimming_axis(2),swimming_axis(3),'linewidth',2);
% %         xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
% %         title(num2str(i / Nt));
% %         view(90,0)
% % %         drawnow
% %         y(7) = mod(y(7) + 2*pi, 2*pi);
% %         temp = cos( acos( sum( swimming_axis .* swimming_axis0 ) ) );
% %         if temp > pi
% %             temp = 2*pi - temp; %angle should
% %         end
%         cos_theta_temp(i) =   sum( swimming_axis .* swimming_axis0 );
%     end
% %     stopa
%     cos_theta(:,j) = cos_theta_temp;
%     if do_progmon
%      ppm.increment(j);
%     end
% end
% if do_progmon
% delete(ppm);
% end


%
%%
mean_cos_theta = mean(cos_theta,2);
[fitobject,gof] = fit(T(1:end)',mean_cos_theta(1:end),'exp1','Lower',[1 -Inf],'Upper',[1 Inf]);
% [fitobject,gof] = fit(T(1:end)',mean_cos_theta(1:end),'exp1','Lower',[1 -0.6135],'Upper',[1 -0.6135]);
figure(66);
plot(T,mean_cos_theta,'o-'); grid on;
hold on;
% plot(T,fitobject.a*exp(fitobject.b*T),'k--','linewidth',2);
pp = csaps(T,mean_cos_theta);
plot(T,fnval(pp,T),'k--','linewidth',2);
hold off

% tau = -1/fitobject.b

tau = fzero( @(t) fnval(pp,t) - 1/exp(1) , 5)
toc
assert(tau < 0.95*T(end),'tau > T(end)');
return
%% translation
X = [repmat(reshape(X0,[1 1 3]),1,Np,1); repmat(reshape(X0,[1 1 3]),Nt,Np,1) + cumsum(dx,1)];

RMSD = squeeze(  sqrt( mean( (X - repmat(reshape(X0,[1 1 3]),Nt+1,Np,1) ).^2 , 2)    )  );

%%

figure(48)
plot(T,RMSD,'o-'); grid on

[fitobject,gof] = fit(T(2:end)',RMSD(2:end),'power1');

hold on;
plot(T,fitobject.a*T.^fitobject.b,'k--','linewidth',2);
hold off

fitobject.b

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
