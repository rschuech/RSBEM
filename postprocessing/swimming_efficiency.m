%
%check for fucked mesh:
%AR1 4
%AR2 0.8

% folder = 'E:\Hull\dumps5\ellipsoid opt\';
%
% folder = 'C:\shatlab\';
plot_datapts = true;
npts = 300;
shape = 'capsule';
%shape = 'ellipsoid';
shape = 'curved_rod';

%clear torque_effs power_effs AR1s AR2s
eff.(shape).torque_effs = [];
eff.(shape).power_effs = [];
eff.(shape).AR1 = [];  eff.(shape).AR2 = [];  eff.(shape).amp = [];  eff.(shape).lambda = [];   eff.(shape).nlambda = [];
eff.(shape).speeds = [];  eff.(shape).powers = [];  eff.(shape).torques = [];  eff.(shape).freqs = [];

AR1_discard = [3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75 7.25 ];
AR2_discard = [0.12 0.14 0.16 0.18 0.22 0.24 0.26 0.28 0.32 0.34 0.36 0.38 ...
               0.42 0.44 0.46 0.48 0.52 0.54 0.56 0.58 0.62 0.64 0.66 0.68 ...
               0.72 0.74 0.76 0.78 0.82 0.84   ];

%folders = {'E:\Hull\dumps5\ellipsoid opt\','C:\shatlab\'};

folders = {'E:\Hull\minidumps\'};  %where the minidump files are with timestepping results
torque_hardcoded = 1E-6;  %very badly hardcoded for now but should eventually read in from interp dump files, or copy from interp dump files into minidump files


filecount = 0;
for fold = 1:length(folders)
    
    folder = folders{fold};
    Files = dir([folder,shape,'*converged_speed.mat']);
    filecount = filecount + length(Files);
end

ftot = 0;
for fold = 1:length(folders)
    
    folder = folders{fold};
    
    Files = dir([folder,shape,'*converged_speed.mat']);
    Files = {Files.name};
    
    
    
    for ff = 1:length(Files)
        ftot = ftot + 1;
        
       % doneness = [ftot   ftot/filecount];
       if rem(ftot, 200) == 0
        disp(['Iter ',num2str(ftot),'     ',num2str(ftot/filecount*100),' % done']);
       end
        
        File = Files{ff};
        
        if strcmp(File,'curved_rod_AR1_10.5_AR2_0.4_tail_radius_0.031018_amp_0.75321_lambda_3.245_nlambda_1.9381_motorBC_torque_converged_speed.mat')
            pause
        end
        
        AR1 = str2double(File(strfind(File,'AR1_')+4 : strfind(File,'_AR2')-1));
        AR2 = str2double(File(strfind(File,'AR2_')+4 : strfind(File,'_tail')-1));
        amp = str2double(File(strfind(File,'amp_')+4 : strfind(File,'_lambda')-1));
        lambda = str2double(File(strfind(File,'lambda_')+7 : strfind(File,'_nlambda')-1));
        nlambda = str2double(File(strfind(File,'nlambda_')+8 : strfind(File,'_motorBC')-1));
        
  if ismember(AR1, AR1_discard) || ismember(AR2, AR2_discard)
      continue
  end
        
        %         if ~isempty(eff.(shape).ARs) && ismember([AR AR2],[eff.(shape).ARs' eff.(shape).AR2'],'rows')
        %             continue
        %         end
        
        %       if strcmp(shape,'capsule') && AR2 ~= 1
        %           continue
        %       end
        
        %         if ~isempty(findstr(File,'fast'))  % fast method
        %                 load([folder,Files{ff}],'c','parameters','solf','n_col','ders','thetas','sol','Geom_B');
        %
        %                 %continue
        %         else
        %                 load([folder,Files{ff}],'c','parameters','sol','solf','n_col');
        %         end
        load([folder,Files{ff}]);
        
        
        eff.(shape).AR1(end+1) = AR1;
        eff.(shape).AR2(end+1) = AR2;
        eff.(shape).amp(end+1) = amp;
        eff.(shape).lambda(end+1) = lambda;
        eff.(shape).nlambda(end+1) = nlambda;
        
        % File
        
        
        %         if  abs(eff.(shape).ARs(end) - parameters(c).body.aspect_ratio) / parameters(c).body.aspect_ratio >= 0.01  ||  abs(eff.(shape).AR2(end) - parameters(c).body.aspect_ratio2) / parameters(c).body.aspect_ratio2 >= 0.01
        %
        %             %  Files{ff}
        %             actual_ARs =      [parameters(c).body.aspect_ratio  parameters(c).body.aspect_ratio2]
        %             filename_ARs = [AR AR2]
        %
        %
        %             eff.(shape).ARs(end) = parameters(c).body.aspect_ratio;
        %             eff.(shape).AR2(end) = parameters(c).body.aspect_ratio2;
        %
        %
        %         end
        
        %only works for constant torque condition right now
        
        avg_freq = timestepping_convergence.avg_omega;
        avg_torque = torque_hardcoded;
        avg_power = avg_freq * avg_torque;  %one of these should be constant and the other varying unless I eventually implement a constant power motor condition
        eff_speed = timestepping_convergence.avg_speed;
        
        mu = 1E-9; %kg / (s micron)
        % effspeed = 18.33;
        % effspeed = 36.51;
        % motor_torque = 1E-6;  %microN micron = kg micron^2 / s^2
        V = 1 ;  %all shapes will have this volume = 1 um^3
        sphererad = (V*3/4/pi)^(1/3); %only used to get tail radius according to Shum et al
        
        
        torque_eff = 8*pi*mu*sphererad^2*eff_speed / avg_torque;
        
        
        % avgpower = motor_torque * 619;
        % avgpower = motor_torque * 14993;
        
        power_eff = 6*pi*mu*sphererad*eff_speed^2 / avg_power;
        %1E3 * power_eff
        
        eff.(shape).torque_effs(end+1) = torque_eff;
        eff.(shape).power_effs(end+1) = power_eff;
        
        eff.(shape).speeds(end+1) = eff_speed;
        eff.(shape).powers(end+1) = avg_power;
        eff.(shape).torques(end+1) = avg_torque;
        eff.(shape).freqs(end+1) = avg_freq;
        
        %                 if strcmp(shape,'curved_rod')
        %                 eff.(shape).AR22s(end+1) = (Geom_B.arclength - 2*Geom.radius) / (Geom_B.radius_curv * 2 * pi);
        %                 end
        
    end
    
end


[~,inds] = sort(eff.(shape).power_effs,'descend');
fields = fieldnames(eff.(shape));
for f = 1:length(fields)
    eff.(shape).(fields{f}) = eff.(shape).(fields{f})(inds);
end
stopooo
max_power_eff = 0.1;  %anything over this assumed to be wrong
inds = find( isnan(eff.(shape).power_effs) | eff.(shape).power_effs > max_power_eff);
for f = 1:length(fields)
    eff.(shape).(fields{f})(inds) = [];
end
%%
AR1s = unique(eff.curved_rod.AR1);
AR2s = unique(eff.curved_rod.AR2);
amps = unique(eff.curved_rod.amp);  amps(2,:) = amps / amps(2);
lambdas = unique(eff.curved_rod.lambda);  lambdas(2,:) = lambdas / lambdas(2);
nlambdas = unique(eff.curved_rod.nlambda);  nlambdas(2,:) = nlambdas / nlambdas(3);
%%
stopafra

% bad_inds = 1;
% fields = fieldnames(eff.curved_rod);
% for f = 1:length(fields)
%     field = fields{f};
%     eff.curved_rod.(field) = eff.curved_rod.(field)(setdiff(1:length(eff.curved_rod.(field)),bad_inds));
% end
%%


%%
  fontsize = 22;
  
        [X,Y] = meshgrid(linspace(min(eff.(shape).AR1),max(eff.(shape).AR1),npts),linspace(min(eff.(shape).AR2),max(eff.(shape).AR2),npts));
        
        amps = unique(eff.(shape).amp);  lambdas = unique(eff.(shape).lambda);  nlambdas = unique(eff.(shape).nlambda);
        
        
        amp = amps(end);    lambda = lambdas(end);    nlambda = nlambdas(2);
       % amp = amps(3);    lambda = lambdas(1);    nlambda = nlambdas(5);
        
        F = scatteredInterpolant(eff.(shape).AR1(eff.(shape).amp == amp & eff.(shape).lambda == lambda & eff.(shape).nlambda == nlambda)',...
            eff.(shape).AR2(eff.(shape).amp == amp & eff.(shape).lambda == lambda & eff.(shape).nlambda == nlambda)',...
            eff.(shape).power_effs(eff.(shape).amp == amp & eff.(shape).lambda == lambda & eff.(shape).nlambda == nlambda)' / max(eff.(shape).power_effs),...
            'linear','none');
        %F = scatteredInterpolant(ARs_mirrored',AR2s_mirrored',torque_effs_mirrored'/max(torque_effs_mirrored),'natural','none');
        
        Z = F(X,Y);
        
        %     figure(456)
        %     surf(X,Y,Z)
        %     xlabel('AR')
        %     ylabel('AR2')
        %     zlabel('power eff')
        %     hold on
        %     d = plot3(ARs_mirrored,AR2s_mirrored,power_effs_mirrored / max(power_effs_mirrored),'o','markerfacecolor','k');
        %     %d = plot3(ARs_mirrored,AR2s_mirrored,torque_effs_mirrored / max(torque_effs_mirrored),'o','markerfacecolor','k');
        %
figure(766)
        pcolor(X,Y,Z)
       % surf(X,Y,Z)
        shading interp
        xlabel('AR_1','fontsize',fontsize)
        ylabel('AR_2','fontsize',fontsize)
        %zlabel('power eff')
        hold on
        con = contour(X,Y,Z,[0.9 0.95],'k-');
        clabel(con);
        if plot_datapts
            d = plot(eff.(shape).AR1(eff.(shape).amp == amp & eff.(shape).lambda == lambda & eff.(shape).nlambda == nlambda),...
                eff.(shape).AR2(eff.(shape).amp == amp & eff.(shape).lambda == lambda & eff.(shape).nlambda == nlambda),...
                '.','markersize',8,'markerfacecolor','none','markeredgecolor','k');
        end
     %   plot([0 20],[0 20],'k--');
        hold off
        %axis equal
        set(gca,'fontsize',fontsize*.8);
        grid on
        cbar = colorbar;
        set(cbar,'fontsize',fontsize*.8);
        %xlim([0 10]);
        ylim([0 1]);
        caxis([0 1]);
        %cbl = cblabel('swimming efficiency / best efficiency','fontsize',fontsize);
        %%



%%


if strcmp(shape,'curved_rod')  %copy circular capsule results as straight AR2 = 0 cases
    
    eff.(shape).ARs = [eff.(shape).ARs  eff.capsule.ARs(eff.capsule.AR2s == 1)];
    eff.(shape).AR2 = [eff.(shape).AR2  zeros(1,sum(eff.capsule.AR2s == 1))];
    eff.(shape).power_effs = [eff.(shape).power_effs  eff.capsule.power_effs(eff.capsule.AR2s == 1)];
    eff.(shape).torque_effs = [eff.(shape).torque_effs  eff.capsule.torque_effs(eff.capsule.AR2s == 1)];
    eff.(shape).speeds = [eff.(shape).speeds  eff.capsule.speeds(eff.capsule.AR2s == 1)];
    eff.(shape).powers = [eff.(shape).powers  eff.capsule.powers(eff.capsule.AR2s == 1)];
    eff.(shape).torques = [eff.(shape).torques  eff.capsule.torques(eff.capsule.AR2s == 1)];
    
    
    [~,inds] = unique([eff.(shape).AR1' eff.(shape).AR2' ],'rows');
    
    eff.(shape).ARs = eff.(shape).ARs(inds);
    eff.(shape).AR2 = eff.(shape).AR2(inds);
    eff.(shape).power_effs = eff.(shape).power_effs(inds);
    eff.(shape).torque_effs = eff.(shape).torque_effs(inds);
    eff.(shape).speeds = eff.(shape).speeds(inds);
    eff.(shape).powers = eff.(shape).powers(inds);
    eff.(shape).torques = eff.(shape).torques(inds);
    
end



%AR2s_back = AR2s;
if ~strcmp(shape,'curved_rod')
    eff.(shape).AR22s = eff.(shape).ARs .* eff.(shape).AR2;
    
    %mirror AR2s around 1 since AR2 = 4 is same as AR2 = 1/4, etc
    %%
    ARs_mirrored = [eff.(shape).ARs eff.(shape).AR22s];
    AR2s_mirrored = [eff.(shape).AR22s eff.(shape).ARs];
    power_effs_mirrored = [eff.(shape).power_effs eff.(shape).power_effs];
    torque_effs_mirrored = [eff.(shape).torque_effs eff.(shape).torque_effs];
    [~,inds] = unique([ARs_mirrored' AR2s_mirrored' power_effs_mirrored'],'rows');
    
    eff.(shape).ARs_mirrored = ARs_mirrored(inds);
    eff.(shape).AR2s_mirrored = AR2s_mirrored(inds);
    eff.(shape).power_effs_mirrored = power_effs_mirrored(inds);
    eff.(shape).torque_effs_mirrored = torque_effs_mirrored(inds);
end

switch shape
    case {'capsule','ellipsoid'}
        
        %% capsule or ellipsoid
        
        fontsize = 22;
        [X,Y] = meshgrid(linspace(min(eff.(shape).ARs_mirrored),max(eff.(shape).ARs_mirrored),npts),linspace(min(eff.(shape).AR2s_mirrored),max(eff.(shape).AR2s_mirrored),npts));
        F = scatteredInterpolant(eff.(shape).ARs_mirrored',eff.(shape).AR2s_mirrored',eff.(shape).power_effs_mirrored'/max(eff.(shape).power_effs_mirrored),'natural','none');
        %F = scatteredInterpolant(ARs_mirrored',AR2s_mirrored',torque_effs_mirrored'/max(torque_effs_mirrored),'natural','none');
        
        Z = F(X,Y);
        
        %     figure(456)
        %     surf(X,Y,Z)
        %     xlabel('AR')
        %     ylabel('AR2')
        %     zlabel('power eff')
        %     hold on
        %     d = plot3(ARs_mirrored,AR2s_mirrored,power_effs_mirrored / max(power_effs_mirrored),'o','markerfacecolor','k');
        %     %d = plot3(ARs_mirrored,AR2s_mirrored,torque_effs_mirrored / max(torque_effs_mirrored),'o','markerfacecolor','k');
        %
        %     hold off
        switch shape
            case 'ellipsoid'
                figure(457)
            case 'capsule'
                figure(458)
        end
        pcolor(X,Y,Z)
        shading interp
        xlabel('AR_1','fontsize',fontsize)
        ylabel('AR_2','fontsize',fontsize)
        %zlabel('power eff')
        hold on
        con = contour(X,Y,Z,0.9,'w-');
        if plot_datapts
            d = plot(eff.(shape).ARs_mirrored,eff.(shape).AR2s_mirrored,'.','markersize',4,'markerfacecolor','none','markeredgecolor','k');
        end
        plot([0 20],[0 20],'k--');
        hold off
        axis equal
        set(gca,'fontsize',fontsize*.8);
        grid on
        cbar = colorbar;
        set(cbar,'fontsize',fontsize*.8);
        xlim([0 10]);
        ylim([0 10]);
        caxis([0.25 1]);
        cbl = cblabel('swimming efficiency / best efficiency','fontsize',fontsize);
        %%
    case 'curved_rodfuck'
        %% curved rod
        
        
        [X,Y] = meshgrid(linspace(min(eff.curved_rod.ARs),max(eff.curved_rod.ARs),npts),linspace(min(eff.curved_rod.AR2s),max(eff.curved_rod.AR2s),npts));
        F = scatteredInterpolant(eff.curved_rod.ARs',eff.curved_rod.AR2s',eff.curved_rod.power_effs'/max(eff.curved_rod.power_effs),'natural','none');
        %F = scatteredInterpolant(ARs_mirrored',AR2s_mirrored',torque_effs_mirrored'/max(torque_effs_mirrored),'natural','none');
        
        Z = F(X,Y);
        
        figure(459)
        pcolor(X,Y,Z)
        shading interp
        xlabel('AR')
        ylabel('AR2')
        hold on
        d = plot(eff.curved_rod.ARs,eff.curved_rod.AR2s,'.','markersize',5,'markerfacecolor','none','markeredgecolor','k');
        hold off
        %axis equal
        grid on
        colorbar
        xlim([1 Inf])
        ylim([0 1])
        caxis([0.5 1]);
        
end

%%

[eff.ellipsoid.ARs' eff.ellipsoid.AR22s' eff.ellipsoid.AR2s' eff.ellipsoid.power_effs' eff.ellipsoid.speeds' eff.ellipsoid.torques' eff.ellipsoid.powers']

[eff.ellipsoid.ARs(eff.ellipsoid.ARs == eff.ellipsoid.AR22s)' eff.ellipsoid.AR22s(eff.ellipsoid.ARs == eff.ellipsoid.AR22s)' eff.ellipsoid.power_effs(eff.ellipsoid.ARs == eff.ellipsoid.AR22s)' eff.ellipsoid.speeds(eff.ellipsoid.ARs == eff.ellipsoid.AR22s)' eff.ellipsoid.torques(eff.ellipsoid.ARs == eff.ellipsoid.AR22s)' eff.ellipsoid.powers(eff.ellipsoid.ARs == eff.ellipsoid.AR22s)']
figure(700); plot(eff.ellipsoid.ARs(eff.ellipsoid.ARs == eff.ellipsoid.AR22s)',  eff.ellipsoid.power_effs(eff.ellipsoid.ARs == eff.ellipsoid.AR22s)','o'); grid


[eff.capsule.ARs' eff.capsule.AR22s' eff.capsule.AR2s' eff.capsule.power_effs' eff.capsule.speeds' eff.capsule.torques' eff.capsule.powers']
[eff.capsule.ARs(eff.capsule.ARs == eff.capsule.AR22s)' eff.capsule.AR22s(eff.capsule.ARs == eff.capsule.AR22s)' eff.capsule.AR2s(eff.capsule.ARs == eff.capsule.AR22s)' eff.capsule.power_effs(eff.capsule.ARs == eff.capsule.AR22s)' eff.capsule.speeds(eff.capsule.ARs == eff.capsule.AR22s)' eff.capsule.torques(eff.capsule.ARs == eff.capsule.AR22s)' eff.capsule.powers(eff.capsule.ARs == eff.capsule.AR22s)']
figure(701); plot(eff.capsule.ARs(eff.capsule.ARs == eff.capsule.AR22s)',  eff.capsule.power_effs(eff.capsule.ARs == eff.capsule.AR22s)' /eff.capsule.power_effs(eff.capsule.ARs == 1 & eff.capsule.AR22s == 1) ,'o'); grid


[eff.curved_rod.ARs' eff.curved_rod.AR2s' eff.curved_rod.power_effs' eff.curved_rod.speeds' eff.curved_rod.torques' eff.curved_rod.powers']

AR1 = 8;
[eff.curved_rod.ARs(eff.curved_rod.ARs == AR1)' eff.curved_rod.AR2s(eff.curved_rod.ARs == AR1)' eff.curved_rod.power_effs(eff.curved_rod.ARs == AR1)' eff.curved_rod.speeds(eff.curved_rod.ARs == AR1)' eff.curved_rod.torques(eff.curved_rod.ARs == AR1)' eff.curved_rod.powers(eff.curved_rod.ARs == AR1)']
figure(460)
plot(eff.curved_rod.AR2s(eff.curved_rod.ARs == AR1),eff.curved_rod.power_effs(eff.curved_rod.ARs == AR1),'o'); grid on;
%%
if do_ratio
    param_limits.capsule = [min(eff.capsule.ARs_mirrored) max(eff.capsule.ARs_mirrored); min(eff.capsule.AR2s_mirrored) max(eff.capsule.AR2s_mirrored)];
    
    param_limits.ellipsoid = [min(eff.ellipsoid.ARs_mirrored) max(eff.ellipsoid.ARs_mirrored); min(eff.ellipsoid.AR2s_mirrored) max(eff.ellipsoid.AR2s_mirrored)];
    
    minlims = [max(param_limits.capsule(:,1),param_limits.ellipsoid(:,1)) min(param_limits.capsule(:,2),param_limits.ellipsoid(:,2)); ];
    
    [X,Y] = meshgrid(linspace(minlims(1,1),minlims(1,2),npts),linspace(minlims(2,1),minlims(2,2),npts));
    F = scatteredInterpolant(eff.capsule.ARs_mirrored',eff.capsule.AR2s_mirrored',eff.capsule.power_effs_mirrored','natural','none');
    Z_capsule = F(X,Y);
    F = scatteredInterpolant(eff.ellipsoid.ARs_mirrored',eff.ellipsoid.AR2s_mirrored',eff.ellipsoid.power_effs_mirrored','natural','none');
    Z_ellipsoid = F(X,Y);
    
    ratio = Z_capsule ./ Z_ellipsoid;
    
    figure(394)
    pcolor(X,Y,ratio)
    shading interp
    xlabel('AR')
    ylabel('AR2')
    %zlabel('power eff')
    hold on
    %  d = plot(eff.(shape).ARs_mirrored,eff.(shape).AR2s_mirrored,'.','markersize',5,'markerfacecolor','none','markeredgecolor','k');
    plot([0 25],[0 25],'k--');
    hold off
    axis equal
    grid on
    colorbar
    
end




%%
fontsize = 22;
% colors = ['kkk'];
color = [repmat(0.4,1,3)];
markers = ['^so'];
% maxeff = max([eff.capsule.power_effs eff.ellipsoid.power_effs eff.curved_rod.power_effs]);
%  maxeff = max([eff.capsule.power_effs  ]);
i = 0;
for AR1 = [3 5 8]
    i = i+1;
    [eff.curved_rod.ARs(eff.curved_rod.ARs == AR1)' eff.curved_rod.AR2s(eff.curved_rod.ARs == AR1)' eff.curved_rod.power_effs(eff.curved_rod.ARs == AR1)' eff.curved_rod.speeds(eff.curved_rod.ARs == AR1)' eff.curved_rod.torques(eff.curved_rod.ARs == AR1)' eff.curved_rod.powers(eff.curved_rod.ARs == AR1)']
    figure(460)
    plot(eff.curved_rod.AR2s(eff.curved_rod.ARs == AR1),eff.curved_rod.power_effs(eff.curved_rod.ARs == AR1)/maxeff,markers(i),'markerfacecolor',color,'markeredgecolor',color);
    
    % plot(eff.curved_rod.AR2s(eff.curved_rod.ARs == AR),eff.curved_rod.speeds(eff.curved_rod.ARs == AR).^2 ./ eff.curved_rod.powers(eff.curved_rod.ARs == AR),markers(i),'markerfacecolor',color,'markeredgecolor',color);
    
    
    hold on
end
grid on;

hold off
xlabel('AR_2','fontsize',fontsize)
ylabel('swimming efficiency / best efficiency','fontsize',fontsize)
set(gca,'fontsize',fontsize*.8);
lh = legend('AR_1 = 3','AR_1 = 5','AR_1 = 8');
set(lh,'fontsize',fontsize*0.7);
ylim([0.65 1]);

hold on
plot(0,eff.capsule.power_effs(find(eff.capsule.ARs == 3 & eff.capsule.AR22s == 3))/maxeff,markers(1),'markerfacecolor',color,'markeredgecolor',color);
plot(0,eff.capsule.power_effs(find(eff.capsule.ARs == 5 & eff.capsule.AR22s == 5))/maxeff,markers(2),'markerfacecolor',color,'markeredgecolor',color);
plot(0,eff.capsule.power_effs(find(eff.capsule.ARs == 8 & eff.capsule.AR22s == 8))/maxeff,markers(3),'markerfacecolor',color,'markeredgecolor',color);
hold off

%%


print(gcf,'capsule efficiency.png','-dpng','-r2000');


% 0.25          2.5            8            6
set(gcf,'paperPosition',[    0.25          2.5            9*1.2            5*1.2]);
print(gcf,'curved rod efficiency.png','-dpng','-r2000');

print(gcf,'ellipsoid efficiency.png','-dpng','-r2000');



%%
npts = 300;
%npts = 20;

var1 = eff.curved_rod.AR1;   var2 = eff.curved_rod.AR2;
%var1 = eff.curved_rod.lambdas;  var2 = eff.curved_rod.nlambdas;
%var1 = eff.curved_rod.nlambdas;  var2 = eff.curved_rod.amps;
X_Y = [[var1]; [var2] ]';  %two vars to use for X and Y axes

[X_Y_unq,inds,inds2] = unique(X_Y,'rows');
best_effs = NaN(size(X_Y_unq,1),1);  
best_inds = best_effs; n_pts = best_effs;  best_speeds = best_effs;  best_freqs = best_effs;
best_amps = best_effs;  best_lambdas = best_effs;  best_nlambdas = best_effs;
n_amps = best_effs;  n_lambdas = best_effs;  n_nlambdas = best_effs;

for i = 1:size(X_Y_unq,1) %each unique body shape
    x_y = X_Y_unq(i,:);
    inds_temp = find(X_Y(:,1) == x_y(1) & X_Y(:,2) == x_y(2)); %all the runs with this body shape
    effs = eff.curved_rod.power_effs(inds_temp);
    best_ind = inds_temp( effs == max(effs));  %index of best run within this subset
    best_effs(i) = eff.curved_rod.power_effs(best_ind);
    best_inds(i) = best_ind;
    n_pts(i) = length(inds_temp);
    best_speeds(i) = eff.curved_rod.speeds(best_ind);
    best_freqs(i) = eff.curved_rod.freqs(best_ind);
    best_amps(i) = eff.curved_rod.amp(best_ind);
    best_lambdas(i) = eff.curved_rod.lambda(best_ind);
    best_nlambdas(i) = eff.curved_rod.nlambda(best_ind);
    n_amps(i) = length(unique(eff.curved_rod.amp(inds_temp)));
    n_lambdas(i) = length(unique(eff.curved_rod.lambda(inds_temp)));
    n_nlambdas(i) = length(unique(eff.curved_rod.nlambda(inds_temp)));
    
end




  [X,Y] = meshgrid(linspace(min(var1),max(var1),npts),linspace(min(var2),max(var2),npts));
        dep_var = best_effs / max(eff.(shape).power_effs);
     %    dep_var =  best_speeds;
%          dep_var = best_freqs;
        %  dep_var = best_nlambdas;
%            dep_var = best_amps;
%           dep_var = best_nlambdas;
 % dep_var = n_pts;
% dep_var = best_speeds ./ best_freqs;
%         

%suspected region of global optimum:
% AR1  3.5 - 7
% AR2  0.45 - 0.75
% amp > 0.5025
% lambda > 3.6291
% nlambda > 1.1175 &  < 1.6763

%suspected region of straight optimum:
% AR1  1 - 3
% AR2  0
% amp > 0.3015 & < 0.5025
% lambda > 2.1774 & < 3.6291
% nlambda > 1.1175 &  < 1.49


      subset = n_pts >= 1;
      method = 'natural';
       method = 'nearest';
        F = scatteredInterpolant(X_Y_unq(subset,1), X_Y_unq(subset,2), dep_var(subset),method,'none');
        %F = scatteredInterpolant(ARs_mirrored',AR2s_mirrored',torque_effs_mirrored'/max(torque_effs_mirrored),'natural','none');
        
        Z = F(X,Y);
        figure(778)
        %     figure(456)
        %     surf(X,Y,Z)
        %     xlabel('AR')
        %     ylabel('AR2')
        %     zlabel('power eff')
        %     hold on
        %     d = plot3(ARs_mirrored,AR2s_mirrored,power_effs_mirrored / max(power_effs_mirrored),'o','markerfacecolor','k');
        %     %d = plot3(ARs_mirrored,AR2s_mirrored,torque_effs_mirrored / max(torque_effs_mirrored),'o','markerfacecolor','k');
        %
fontsize = 22;
        pcolor(X,Y,Z)
       % surf(X,Y,Z)
        shading interp
        xlabel('AR_1','fontsize',fontsize)
        ylabel('AR_2','fontsize',fontsize)
%         
%             xlabel('tail amplitude (\mum)','fontsize',fontsize)
%         xlabel('tail wavelength (\mum)','fontsize',fontsize)
%           ylabel('# tail wavelengths','fontsize',fontsize)
        %zlabel('power eff')
        hold on
       % con = contour(X,Y,Z,[0.98 0.98],'k-');
       % clabel(con);
        if plot_datapts
           try
            d1 = plot(X_Y_unq(subset & isnan([next_iter.amp])',1),X_Y_unq(subset & isnan([next_iter.amp])',2), '.','markersize',8,'markerfacecolor','none','markeredgecolor','k');
            d2 = plot(X_Y_unq(subset & ~isnan([next_iter.amp])',1),X_Y_unq(subset & ~isnan([next_iter.amp])',2), 'o','markersize',5,'markerfacecolor','none','markeredgecolor','k');
           catch
               d = plot(X_Y_unq(subset,1),X_Y_unq(subset,2), '.','markersize',8,'markerfacecolor','none','markeredgecolor','k');
           end
        
        end
     %   plot([0 20],[0 20],'k--');
        hold off
        %axis equal
        set(gca,'fontsize',fontsize*.8);
        grid on
        cbar = colorbar;
        set(cbar,'fontsize',fontsize*.8);
        %xlim([0 10]);
       % ylim([0 1]);
     %  caxis([0.6 1]);
       % cbl = cblabel('swimming efficiency / best efficiency','fontsize',fontsize);
        
        l = legend([d1 d2],'tail crudely optimized','tail not optimized');
        set(l,'position',[0.58521      0.90343      0.23778     0.067114]);
        
         cbl = cblabel('swimming speed (\mum/s)','fontsize',fontsize);
        % cbl = cblabel('motor frequency (rev/s)','fontsize',fontsize);
       %  cbl = cblabel('# tail wavelengths','fontsize',fontsize);
          %cbl = cblabel('tail amplitude (\mum)','fontsize',fontsize);
     %    cbl = cblabel('tail wavelength (\mum)','fontsize',fontsize);
     %  cbl = cblabel('# tails tested','fontsize',fontsize);
%         best_line.AR1 = X(1,:);  %all AR1s
%         best_line.AR2 = NaN(size(X(1,:)));  %best corresponding AR2s to be filled in below
%         for i = 1:size(X,2)
%             [~, max_eff_ind] = max(Z(:,i));
%             best_line.AR2(i) = Y(max_eff_ind,1);
%         end
%         hold on
%         pl = plot(best_line.AR1,best_line.AR2,'k-');
%         hold off
            
        %%
        
   varlist = {'AR1s','AR2s','amps','lambdas','nlambdas'};     
        % look at best choices of two vars, given each unique combo of all
        % other vars
        variable1 = 'AR1s';  variable2 = 'AR2s';
var1 = eff.curved_rod.(variable1);   var2 = eff.curved_rod.(variable2);  
%var1 = eff.curved_rod.amps;  var2 = eff.curved_rod.lambdas;
others = setdiff(varlist,{variable1, variable2});

othervars = [];
for i = 1:length(others)
othervars = [othervars  eff.curved_rod.(others{i})'];  %two vars to use for X and Y axes
end

[othervars_unq,inds,inds2] = unique(othervars,'rows');
% best.effs = NaN(size(othervars_unq,1),1);  
fields = fieldnames(eff.curved_rod);
for f = 1:length(fields)
    field = fields{f};
    best.(field) = NaN(size(othervars_unq,1),1);  
end
best.ind = NaN(size(othervars_unq,1),1);  
best.npts = NaN(size(othervars_unq,1),1);  
% best_var1s = best_effs;
% best_var2s = best_effs;

for i = 1:size(othervars_unq,1) %each unique combo of other vars
    othervar = othervars_unq(i,:);
    logical_temp = true(size(othervars,1),1);
    for j = 1:size(othervars,2)  %each var type, usually 3
        logical_temp = logical_temp & othervars(:,j) == othervar(j);  %incremental AND aggregation - if any is false, it goes false
    end
    inds_temp = find(logical_temp);
    effs = eff.curved_rod.power_effs(inds_temp);
    best_ind = inds_temp( effs == max(effs));  %index of best run within this subset
   % best_effs(i) = eff.curved_rod.power_effs(best_ind);
   for f = 1:length(fields)
       best.(fields{f})(i) = eff.curved_rod.(fields{f})(best_ind);
   end
    best.ind(i) = best_ind;
    best.npts(i) = length(inds_temp);
    
end

%%
figure;  %plot(best.AR1s,best.AR2s,'o','markerfacecolor','k'); grid on
clear text
minpts = 5;
colorlo = 0.65;  colorhi = 1;
[colors,colorinds] = colordata(50,'parula',[colorlo colorhi],best.power_effs(best.npts > minpts) / max(best.power_effs));

x = best.AR1s(best.npts > minpts);  y = best.AR2s(best.npts > minpts);
for i = 1:length(x)
plot(x(i),y(i),'o','markerfacecolor',colors(colorinds(i),:),'markeredgecolor',colors(colorinds(i),:),'markersize',6)
if sum( x == x(i) & y == y(i)) > 1
    text(x(i)+0.01,y(i)+0.01,num2str(sum( x == x(i) & y == y(i))));
end
hold on
end
hold off

colorbar
caxis([colorlo colorhi]);
colormap(colors);
grid on
xlim([1 12]);  ylim([0 1]);
xlabel('AR1'); ylabel('AR2');
title('Best bodies for each unique tail tested')
cblabel('relative efficiency','fontsize',12)
%%
results.power_effs = NaN(length(AR1s),length(AR2s),length(amps(1,:)),length(lambdas(1,:)),length(nlambdas(1,:)));
results.speeds = results.power_effs;
results.freqs = results.speeds;
fields = fieldnames(results);

for i = 1:length(eff.curved_rod.AR1s)
    %     ind = find(AR1s == eff.curved_rod.AR1s(i),AR2s == eff.curved_rod.AR2s(i),amps == eff.curved_rod.amps(i),lambdas == eff.curved_rod.lambdas(i),nlambdas == eff.curved_rod.nlambdas(i));
    %     if length(ind) > 1
    %         stopafra
    %     end
    for f = 1:length(fields)
        results.(fields{f})(AR1s == eff.curved_rod.AR1s(i),AR2s == eff.curved_rod.AR2s(i),amps(1,:) == eff.curved_rod.amps(i),lambdas(1,:) == eff.curved_rod.lambdas(i),nlambdas(1,:) == eff.curved_rod.nlambdas(i)) = eff.curved_rod.(fields{f})(i);
    end
end

%%
min_tails = 29;
for i = 1:size(results.power_effs,1) %AR1
    for j = 1:size(results.power_effs,2) %AR2
        temp = ~isnan(results.power_effs(i,j,:,:,:));
        if sum(temp(:)) >= min_tails
            disp(['AR1 = ',num2str(AR1s(i)),'   AR2 = ',num2str(AR2s(j)),'   has ',num2str(sum(temp(:))),' tails']);
        end
    end
end
%%
colors = distinguishable_colors(20);
AR1 = 4.5;  AR2 = 0.625;
figure(45)
c = 0; clear legends
for i = 1:size(lambdas,2)
    for j = 1:size(nlambdas,2)
        
        if sum(~isnan(squeeze(results.power_effs(AR1s == AR1, AR2s == AR2, :,i,j)))) >= 0
            c = c+1;
            temp = squeeze(results.power_effs(AR1s == AR1, AR2s == AR2, :,i,j));
            x = amps(1,~isnan(temp));  y = temp(~isnan(temp));
            plot(x,y/max(results.power_effs(:)),'-o','markeredgecolor',colors(c,:),'markerfacecolor',colors(c,:),'color',colors(c,:));
            legends{c} = ['\lambda = ',num2str(lambdas(i)),',     # \lambda = ',num2str(nlambdas(j))];
            hold on
            grid on
            
        end
    end
end
hold off
xlabel('tail amplitude (\mum)','fontsize',fontsize);
ylabel('swimming efficiency / best efficiency','fontsize',fontsize);
legend(legends,'fontsize',fontsize-8);
title(['AR_1 = ',num2str(AR1),'       AR_2 = ',num2str(AR2)],'fontsize',fontsize);
        set(gca,'fontsize',fontsize*.8);


%%
%%
sweeps_already_done = [eff.curved_rod.AR1s' eff.curved_rod.AR2s' eff.curved_rod.amps' eff.curved_rod.lambdas' eff.curved_rod.nlambdas'];
