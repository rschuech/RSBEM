tempr = load('C:\Hull\CFD03\dino_code\pape figs\chemotactic SNR vs tail length\Results_lambda.mat');
varied.lambda.tail_parameter = [tempr.Results.lambda];
varied.lambda.Adj_Speed = [tempr.Results.Adj_Speed];
varied.lambda.tau_a = [tempr.Results.tau_a];
temp = [tempr.Results.taxis];  temp = [temp.temporal];  varied.lambda.Chemotactic_SNR = [temp.SN];

tempr = load('C:\Hull\CFD03\dino_code\pape figs\chemotactic SNR vs tail length\Results_nlambda.mat');
varied.nlambda.tail_parameter = [tempr.Results.nlambda];
varied.nlambda.Adj_Speed = [tempr.Results.Adj_Speed];
varied.nlambda.tau_a = [tempr.Results.tau_a];
temp = [tempr.Results.taxis];  temp = [temp.temporal];  varied.nlambda.Chemotactic_SNR = [temp.SN];

[~,ind] = ismember([1 0],Best.temporal.body,'rows');
normalizer.Adj_Speed = Best.temporal.Adj_Speed(ind);
normalizer.tau_a = Best.temporal.tau_a(ind);
normalizer.Chemotactic_SNR = F.temporal.F.metric(standardize([1 0],limits));

clear normalized
var1s = {'lambda','nlambda'};
    var2s = {'Adj_Speed','tau_a','Chemotactic_SNR'};
    
for v1 = 1:length(var1s)
    var1 = var1s{v1};
    for v2 = 1:length(var2s)
        var2 = var2s{v2};
normalized.(var1).(var2) = varied.(var1).(var2) ./ normalizer.(var2);
    end
end

fontsize = 18;
figure(713);  clf;  set(gcf,'Color','w');
yyaxis left
[h_speed] = plot(varied.lambda.tail_parameter,normalized.lambda.Adj_Speed,'b--v',     varied.nlambda.tail_parameter,normalized.nlambda.Adj_Speed,'b--^');
h_speed(1).MarkerFaceColor = 'b';  h_speed(2).MarkerFaceColor = 'b';
ylabel(['\fontsize{',num2str(fontsize),'}{0}$\mathrm{\quad}\overline{U}^{r,n}_{chemo}$'],'Interpreter','latex');
set(gca,'FontSize',fontsize);
yyaxis right
h_tau_SNR = plot(varied.lambda.tail_parameter,normalized.lambda.tau_a, 'r--v' ,    varied.nlambda.tail_parameter,normalized.nlambda.tau_a,'r--^',   ...
     varied.lambda.tail_parameter,normalized.lambda.Chemotactic_SNR,'r-v' ,     varied.nlambda.tail_parameter,normalized.nlambda.Chemotactic_SNR,'r-^'       );
set( h_tau_SNR([3 4]) , 'linewidth', 2);  set(h_tau_SNR,'markerfacecolor','r');
% ylabel({'\tau^{\itn}_{\itchemo}', '\Psi_{\itchemo}'});
% ylabel({['\fontsize{',num2str(fontsize),'}{0}','$\tau^{n}_{chemo}$'], '$\smallskip$' , ['\fontsize{',num2str(fontsize),'}{0}',   '$\Psi_{chemo}$']},'Interpreter','latex');
ylabel(['\fontsize{',num2str(fontsize),'}{0}','$\tau^{n}_{chemo}$' ,'\thinspace \fontsize{11}{0}  \textrm{    or    } \thinspace', '\fontsize{',num2str(fontsize),'}{0}',   '$\Psi_{chemo}$'],'Interpreter','latex');
% xlabel({['\fontsize{',num2str(fontsize),'}','\lambda (\mum)'], ['\fontsize{',num2str(fontsize),'}','\itn_\lambda']},'fontsize',fontsize)
xlabel(['\fontsize{',num2str(fontsize+3),'}','\lambda (\mum)','    or    ','\itn_\lambda'], 'fontsize',fontsize+3);
% set(gca,'FontSize',fontsize);

% export_fig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\SNR_vs_tail_length.png','-r500','-nocrop')

