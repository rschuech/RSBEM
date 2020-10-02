% Coplanar_length = [0:0.5:4];
Coplanar_length = [0:0.5:3];
Normal_length = [0:0.5:1.5];
folder = 'C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs_hair_sweep\';
[Coplanar_Length, Normal_Length] = ndgrid(Coplanar_length,Normal_length);

Speed = NaN(size(Coplanar_Length));  Omega1 = Speed;  Omega2 = Omega1;

for i = 1:numel(Coplanar_Length)
    temp = dir([folder,num2str(Coplanar_Length(i)),'_',num2str(Normal_Length(i)),'_Body-Transverse-Tail*.fig']);
        name = {temp.name};
        if numel(name) > 1
        stopa
        end
       
%     name = [num2str(Coplanar_Length(i)),'_',num2str(Normal_Length(i)),'_Body-Transverse-Tail-Coplanar_Hairs_interpolation.fig'];
    if ~isempty(name) && exist([folder,name{1}],'file')
    uiopen([folder,name{1}],1)
    
    
    clear sols
    chil = get(gcf,'Children');
    for j = 1:6
        sols(1,j) = chil(j).Children.YData;
    end
    
    close
    else
        sols = NaN(1,6);
    end
    
    Speed(i) = sols(6);  Omega1(i) = sols(3);  Omega2(i) = sols(2);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Speed_tailess = NaN(size(Coplanar_Length));  Omega1_tailess = Speed_tailess;  Omega2_tailess = Speed_tailess;
folder = 'C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs_hair_sweep\no tail\';
for i = 1:numel(Coplanar_Length)
    temp = dir([folder,num2str(Coplanar_Length(i)),'_',num2str(Normal_Length(i)),'_Body-Transvers*.fig']);
        name = {temp.name};
        if numel(name) > 1
        stopa
        end
       
%     name = [num2str(Coplanar_Length(i)),'_',num2str(Normal_Length(i)),'_Body-Transverse-Tail-Coplanar_Hairs_interpolation.fig'];
    if ~isempty(name) && exist([folder,name{1}],'file')
    uiopen([folder,name{1}],1)
    
    
    clear sols
    chil = get(gcf,'Children');
    for j = 1:6
        sols(1,j) = chil(j).Children.YData;
    end
    
    close
    else
        sols = NaN(1,6);
    end
    
    Speed_tailess(i) = sols(6);  Omega1_tailess(i) = sols(3);  Omega2_tailess(i) = sols(2);
    
end





%% Speed
method = 'linear';
F = scatteredInterpolant(Coplanar_Length(:), Normal_Length(:), Speed(:),method,'none');
cl = linspace(min(Coplanar_length),max(Coplanar_length),200);  nl = linspace(min(Normal_length),max(Normal_length),200);
[CL,NL] = ndgrid(cl,nl);
Speed_interp = F(CL,NL);

figure(510);  set(gcf,'Position',[317         547        1537         544]);
pcolor(CL,NL,Speed_interp);  shading interp;  colorbar;
xlabel('Coplanar hair length (\mum)');  ylabel('Normal hair length (\mum)');
hold on
  [C,ch] = contour(CL,NL,Speed_interp,round(linspace(45,310,15)),'k--','linewidth',1);
clabel(C,ch,'LabelSpacing',500,'fontsize',12);
plot(Coplanar_Length(~isnan(Speed)),Normal_Length(~isnan(Speed)),'ko','markerfacecolor','k','markersize',3);
hold off
set(gca,'fontsize',14)
cblabel('swimming speed (\mum/s)','fontsize',16)
hold on
 [C_obs,ch_obs] = contour(CL,NL,Speed_interp,[108 108],'k-','linewidth',2);
 obs_str = text(1.3,0.75,'observed speed');
 obs_str.Position = [1.37 0.95 0];  obs_str.Rotation = -85;   obs_str.FontWeight = 'bold';   obs_str.FontSize = 12;
%  clabel(C_obs,ch_obs,'LabelSpacing',800);
 hold off
 
 export_fig(gcf,'C:\Users\rudi\Desktop\RD\pape main results\figures\hair length sweep speed.png','-r400','-transparent')
%% Body rotation
method = 'linear';
F = scatteredInterpolant(Coplanar_Length(:), Normal_Length(:),Omega1(:),method,'none');
cl = linspace(min(Coplanar_length),max(Coplanar_length),200);  nl = linspace(min(Normal_length),max(Normal_length),200);
[CL,NL] = ndgrid(cl,nl);
Omega1_interp = F(CL,NL) * 1/2/pi;

figure(511);  set(gcf,'Position',[317         547        1537         544]);
pcolor(CL,NL,Omega1_interp);  shading interp;  colorbar;
xlabel('Coplanar hair length (\mum)');  ylabel('Normal hair length (\mum)');
hold on
  [C,ch] = contour(CL,NL,Omega1_interp,[-0.5 -0.25 -0.2 -0.1  0.1 0.25 0.4 0.5 ],'k--','linewidth',1);  % roundn(linspace(-0.8,3.5,65),-2)
clabel(C,ch,'LabelSpacing',500,'fontsize',12);
 [C,ch] = contour(CL,NL,Omega1_interp,[0 0],'k--','linewidth',2);  % roundn(linspace(-0.8,3.5,65),-2)
clabel(C,ch,'LabelSpacing',500,'fontsize',12,'fontweight','bold');
plot(Coplanar_Length(~isnan(Speed)),Normal_Length(~isnan(Speed)),'ko','markerfacecolor','k','markersize',3);
hold off
cblabel('rotation rate (rev/s)','fontsize',16)
set(gca,'fontsize',14)
    
export_fig(gcf,'C:\Users\rudi\Desktop\RD\pape main results\figures\hair length sweep rotation.png','-r400','-transparent')

%% Body rotation / Speed
method = 'linear';
F = scatteredInterpolant(Coplanar_Length(:), Normal_Length(:),Omega1(:) ./ Speed(:),method,'none');
cl = linspace(min(Coplanar_length),max(Coplanar_length),200);  nl = linspace(min(Normal_length),max(Normal_length),200);
[CL,NL] = ndgrid(cl,nl);
Ratio_interp = F(CL,NL) * 1/2/pi  * 1E3;

figure(513)
pcolor(CL,NL,Ratio_interp);  shading interp;  colorbar;
xlabel('Coplanar hair length (\mum)');  ylabel('Normal hair length (\mum)');
hold on
  [C,ch] = contour(CL,NL,Ratio_interp,[-4:0.5:2 1.25 1.7 1.8 ],'k--','linewidth',1);  % roundn(linspace(-0.8,3.5,65),-2)
clabel(C,ch,'LabelSpacing',500);
plot(Coplanar_Length(~isnan(Speed)),Normal_Length(~isnan(Speed)),'ko','markerfacecolor','k','markersize',3);
hold off
cblabel('rotation rate / speed  10^3 * (rev / \mum)')
set(gca,'fontsize',14)

%%
method = 'linear';
F = scatteredInterpolant(Coplanar_Length(:), Normal_Length(:), Speed(:) - Speed_tailess(:),method,'none');
cl = linspace(min(Coplanar_length),max(Coplanar_length),200);  nl = linspace(min(Normal_length),max(Normal_length),200);
[CL,NL] = ndgrid(cl,nl);
Speed_diff_interp = F(CL,NL);

figure(512)
pcolor(CL,NL,Speed_diff_interp);  shading interp;  colorbar;
xlabel('Coplanar hair length (\mum)');  ylabel('Normal hair length (\mum)');
hold on
  [C,ch] = contour(CL,NL,Speed_diff_interp,[4:1:14 14.5],'k--','linewidth',1);
clabel(C,ch,'LabelSpacing',500);
plot(Coplanar_Length(~isnan(Speed)),Normal_Length(~isnan(Speed)),'ko','markerfacecolor','k','markersize',3);
hold off
cblabel('(speed with tail) - (speed without tail) (\mum/s)')
set(gca,'fontsize',14)
    