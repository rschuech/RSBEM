



 width = 0.71;
lengths = [0.8051	0.9234	1.0592	1.2148	1.3934	1.5981	1.833	2.1024	2.4114	2.7658	3.1723	3.6385	4.1733	4.7866	5.4901	6.297	7.2224	8.2839	9.5014	10.8978	12.4994	14.3365	16.4435	18.8602	21.632	24.8113	28.4578	32.6402	37.4373	40];
cyl_height = lengths - width;  %total length - 2*radius

AR1 = lengths / width;
V = pi * (width/2)^2 * cyl_height + 4/3*pi*(width/2)^3;

       
       
folder = 'C:\Users\rudi\Desktop\RD\Oscar_dumps\';
clear D fcoeffs D_tailed fcoeffs_tailed


for i = 1:length(AR1)
    i
    name = ['curved_rod_AR1_',num2str(AR1(i),digits),'_AR2_0_forced_dump.mat'];
    name_tailed = ['curved_rod_AR1_',num2str(AR1(i),digits),'_AR2_0_tail_radius_0.012_amp_0.2_lambda_2.34_nlambda_2.564102564102564_forced_dump.mat']; 
    
    temp = load([folder,name]);
    
    D(i) = temp.D;
    fcoeffs(i) = temp.fcoeffs;
    
    if exist([folder,name_tailed],'file')
    temp = load([folder,name_tailed]);
    
    D_tailed(i) = temp.D;
    fcoeffs_tailed(i) = temp.fcoeffs;
    end
    
end
    %%
clear tau tau_ellipsoid tau_tailed
for i = 1:length(AR1)
    tau(i) = 1/(D(i).rotation.diffusivity(2,2) + D(i).rotation.diffusivity(3,3));
    
    if length(D_tailed) >= i
        tau_tailed(i) = 1/(D_tailed(i).rotation.diffusivity(2,2) + D_tailed(i).rotation.diffusivity(3,3));
    end
    
    [fcoeff,frcoeff,tau_a] = ellipsoid_theoretical(lengths(i)/2,width/2,width/2);
    tau_ellipsoid(i) = tau_a;
end
%%
T = 306;  T_0 = 298;
mu = 0.8E-3 / 1E6; %viscosity in kg / (s micron)
mu_0 = 1E-3 / 1E6; %viscosity in kg / (s micron)
correction = 1/(T/T_0*mu_0/mu);
tau_corrected = tau * correction;
tau_ellipsoid_corrected = tau_ellipsoid * correction;
tau_tailed_corrected = tau_tailed * correction;

figure(346)
plot(lengths,tau_corrected,'-ro','markerfacecolor','r','linewidth',2);
hold on
plot(lengths,tau_ellipsoid_corrected,'-bo','markerfacecolor','b','linewidth',2);
plot(lengths(1:length(tau_tailed)),tau_tailed_corrected,'-ks','markerfacecolor','w','linewidth',2)
grid on
hold off
legend('capsule body (numerical)','ellipsoid body (theoretical)','capsule body + tail (numerical)','location','best')
xlabel('body length (\mum)')
ylabel('\tau_{-a} (s)')
set(gca,'fontsize',16)
set(gca,'YScale','log')
% set(gca,'XScale','log');

