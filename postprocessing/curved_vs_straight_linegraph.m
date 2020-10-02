%% find best SF2 for constant SF1 for all metrics


SF1 = linspace(1,10,100);
guesses = unique([0 0.01 0.02 linspace(0,0.94,10)]);
options = optimoptions('patternsearch','meshtolerance',1E-8,'steptolerance',1E-8,'display','none');
clear SF2_best;  SF2_best.SF1 = SF1;

for i = [3:6 ] %1:length(Fs)
    
    
    F_temp = Fs{i};  F_temp.Method = 'natural';  clear SF2_best_temp
    parfor s = 1:length(SF1)
        
        fun = @(SF2) - F_temp( standardize( [SF1(s) SF2] , limits ) );
        best_SF2s = []; fvals = [];
        for g = 1:length(guesses)
            if ~isnan( fun(guesses(g)) )
                [best_SF2s(g) , fvals(g)] = patternsearch(fun,guesses(g),[],[],[],[],0,1,[],options);
            else
                best_SF2s(g) = NaN;  fvals(g) = NaN;
            end
        end
        [~,ind] = min(fvals);
        SF2_best_perf_temp(:,s) = [best_SF2s(ind); -fvals(ind)];
    end
    
    %     [SF1; SF2_best]
    SF2_best.(metrics{i}) = SF2_best_perf_temp;
    
end

%%
% get performances along transect of constant SF2
clear SF2_transect
SF2_transect.SF1 = SF1;
names = {'curved','straight'};
for n = 1:length(names)
    switch names{n}
        case 'curved'
            temp = vertcat(data_filtered.good.median_unweighted);   SF_temp = vertcat(temp.SF);  med = median(SF_temp);
            SF2 = med(2);
            %             avg_AR1_AR2 = mean( vertcat(species_data.SF) );  SF2 = avg_AR1_AR2(2);  %SF2 = 0.4
        case 'straight'
            SF2 = 0;
    end
    
    
    temp = standardize( [SF1'  repmat(SF2,length(SF1),1)] , limits);
    
    for i = [1:6 22]
        SF2_transect.(names{n}).(metrics{i}) = Fs{i}(temp);
    end
    
end

%hardcode known optimal SF2 for certain tasks
for i = [1 2 22] %uptake fore-aft construction
    SF2_best.(metrics{i}) = [zeros(size(SF1));  SF2_transect.straight.(metrics{i})'];
end
%% make line graph of curved vs straight performance
figure(300);  clf;  set(gcf,'color','w');

names = {'tumbling','eff','temporal','uptake','construction'};
% legnames = {'Dispersal','Swimming Efficiency','Chemotaxis (temporal)','Nutrient Uptake','Chemotaxis (spatial)','Construction Ease'};
legnames = {'Tumbling Ease','Swimming Efficiency','Chemotactic SNR','Nutrient Uptake','Construction Ease'};
colors = distinguishable_colors(length(names));
styles = {'-','-','-','--','--','--'};
clear transect_pl
for n = 1:length(names)
    if strcmp('tumbling',names{n})
       [ tpp] = csaps(SF2_transect.SF1 , SF2_transect.curved.(names{n}) ./ SF2_transect.straight.(names{n}), 0.999);
        SF1_refined = linspace(SF2_transect.SF1(2) ,10,200);
        transect_pl(n) = plot(SF1_refined, fnval(tpp,SF1_refined),styles{n});
    else
        
        transect_pl(n) = plot(SF2_transect.SF1 , SF2_transect.curved.(names{n}) ./ SF2_transect.straight.(names{n}),styles{n});
    end
    
    
    transect_pl(n).Color = colors(n,:);  transect_pl(n).LineWidth = 4;
    hold on
end

% names_best = {'eff','Dm','temporal'};  inds = [3 4 5];
% styles_best = {'-','-','-'};
% clear best_pl
% for n = 1:length(names_best)
%     best_pl(n) = plot(SF2_best.SF1 , SF2_best.(names_best{n})(2,:)' ./ SF2_transect.straight.(names_best{n}),styles_best{n});
%     best_pl(n).Color = colors(inds(n),:);  best_pl(n).LineWidth = 3;
%     hold on
% end
% title(['\fontsize{16}black {\color{magenta}magenta ','\color[rgb]{0 .5 .5}teal \color{red}red} black again'],'interpreter','tex')
% title( sprintf('%s%s{%f %f %f}', '\fontsize{16}', '\color[rgb]', RGB_vector), 'interpreter, 'tex')



% legend(transect_pl,legnames,'location','best','interpreter','none','fontsize',12);
legendflex(transect_pl,legnames,'interpreter','none','fontsize',12,'xscale',1.5,'anchor',[1 1],'buffer',[20 -20],'fontsize',18);
hold off
grid on
xlabel('$$ \mathcal{L} $$ (elongation)','interpreter','latex');
% ylabel({'relative performance of','\bf{curved} \rm(SF_2 = 0.15) \bf{vs straight}'});
ylabel({'relative performance of curved vs straight rods'});
set(gca,'fontsize',22);

set(gcf,'Position',[287          84        1355         899]);
set(gcf,'color','w');