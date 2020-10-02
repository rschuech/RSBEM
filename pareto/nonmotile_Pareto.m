A = 1;
center = [0.5 0.6];
width = [1 0.4];

% F is defined in standardized coords
x = linspace(0,1,150);  y = linspace(0,1,150);
[X,Y] = ndgrid(x,y);

C = A * exp( -(  1/2*( (X - center(1))./width(1)   ).^2  +  1/2*( (Y - center(2))./width(2)   ).^2  )  );

figure(349);  pcolor(XP,YP,C); shading interp
F = scatteredInterpolant(X(:), Y(:), C(:),method,'none');

Pareto_data.test1.F = F;

%%

A = 1;
center = [2 -0.5];
width = [0.6 2];

% F is defined in standardized coords
x = linspace(0,1,150);  y = linspace(0,1,150);
[X,Y] = ndgrid(x,y);

C = A * exp( -(  1/2*( (X - center(1))./width(1)   ).^2  +  1/2*( (Y - center(2))./width(2)   ).^2  )  );

figure(350);  pcolor(XP,YP,C); shading interp
F = scatteredInterpolant(X(:), Y(:), C(:),method,'none');

Pareto_data.test2.F = F;


%%


figure(823);  c = 0;   clear legends
for AR1 = [2:2:8]
    c = c+1;
AR1_std = AR1 / (10 - 1); 
AR2_std = linspace(0,1,200); 
transect = Pareto_data.Dbr.F(repmat(AR1_std,size(AR2_std)),AR2_std); 
n = Pareto_data.Dbr.F(AR1_std,0);   % what to normalize by - perf at this AR1 and AR2 = 0
plot((1-0)*AR2_std + 0,(transect - n)/n*100,'-','linewidth',3); hold on;
legends{c} = ['AR_1 = ',num2str(AR1)];
end
hold off;

legend(legends,'location','best')
grid on;

xlabel('AR_2');
ylabel('% improvement over straight');
title('Brownian dispersal')