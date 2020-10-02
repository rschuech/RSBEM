


%%
% min_GOF = [0. 0. 0.];
min_GOF = [0.5 0.0]; %based on [ species individuals] boundaries

nSF1 = 150;  nSF2 = 150;

nSF1 = 50;  nSF2 = 60;



% nSF1 = 100;  nSF2 = 100;
dA = (10 - 1)/(nSF1 - 1) * (1 - 0)/(nSF2 - 1);  %area of a rectangle centered on each interrogation point

[XP,YP] = ndgrid(linspace(1,10,nSF1),linspace(0,1,nSF2));
[standardized_XPYP] = standardize([XP(:) YP(:)], [1 10; 0 1]);
standardized_XP = NaN(size(XP));  standardized_XP(:) = standardized_XPYP(:,1);
standardized_YP = NaN(size(YP));  standardized_YP(:) = standardized_XPYP(:,2);


Pareto_data.interped = NaN(nSF1,nSF2,length(Fs));     Pareto_data.rank =  NaN(nSF1,nSF2,length(Fs));

for f = 1:length(Fs)
    Pareto_data.interped(:,:,f) = Fs{f}(standardized_XP,standardized_YP);
    temp = Pareto_data.interped(:,:,f);
    [~,I] = sort( temp(:) ,'ascend' );
    inds = 1:length(I);
    temp = Pareto_data.rank(:,:,f);
    temp(I) = inds;  %I rearranges inds so that we get the sorted ranks
    
    temp(isnan(Pareto_data.interped(:,:,f))) = NaN;
    Pareto_data.rank(:,:,f) = temp;
end
%%

% goal_inds = setdiff(1:length(Pareto_data.fullnames) , [1 3]);
goal_inds = setdiff(1:length(Pareto_data.fullnames) , []);
% goal_inds = [1:6 18 19];
construction_inds = 7:length(Pareto_data.fullnames);

goal_combos = {};
% loop = length(goal_inds):-1:3;
loop = 4:-1:3;
 loop = [  6 5 4 3 ];

for ngoals = loop
    
    
    temp = combnk(goal_inds,ngoals);
    
    
    for t = 1:size(temp,1)
        
        % don't include more than 2 construction ease variants
        ism = ismember(construction_inds, temp(t,:));
        if sum(ism) > 2
            continue
        end
        
        goal_combos = [goal_combos; {temp(t,:)}];

    end
end

goal_combos = flipud(goal_combos);

%  goal_combos = {[5 6 18]};
% goal_combos = {[5 6 19]};
% 
% 
% goal_combos = {[5 6 7],[5 6 8],[5 6 9],[5 6 10],[5 6 11],[5 6 12],[5 6 13],[5 6 14],[5 6 15],[5 6 16],[5 6 17],[5 6 18],[5 6 19],[5 6 20]};
%%

measurements_only = false;
plot_obs = true;


fontsizes.labels = 21;
fontsizes.axes = 21;
fontsizes.title = 15;
fontsizes.caxis = 26;
fontsizes.legend = 21;

figure(52)
clf; clear pc


if plot_obs
    
    individuals = plot([data.SF1],[data.SF2],'.','markerfacecolor',repmat(0.,1,3),'markeredgecolor',repmat(0.,1,3),'markersize', 5);
    hold on
    temp = vertcat(species_data.mean_unweighted);
    species = plot(temp(:,1),temp(:,2),'o','markerfacecolor','b','markersize',7);
    hold on
%     temp = vertcat(genus_data.mean_unweighted);
%     genera = plot(temp(:,1),temp(:,2),'ro','markerfacecolor','r','markersize',7);
    
%      bound_individuals = plot(Pareto_data.observed.individuals.boundary(:,1),Pareto_data.observed.individuals.boundary(:,2),':','linewidth',1.5,'color',repmat(0.4,1,3));
    bound_species = plot(Pareto_data.observed.species.boundary(:,1),Pareto_data.observed.species.boundary(:,2),'b:','linewidth',1.5);
%     bound_genera = plot(Pareto_data.observed.genera.boundary(:,1),Pareto_data.observed.genera.boundary(:,2),'r-','linewidth',2);
    
    xlabel('SF_1 (elongation)','fontsize',fontsizes.labels)
    ylabel('SF_2 (curvature)','fontsize',fontsizes.labels)
    set(gcf,'Position',[  345         189        1401         729]);
    %         ptemp = patch('Faces',[1 2 3 4],'Vertices',[1 0.9; 1 1; 2 1; 2 0.9],'visible','off');
    hold off
    
end


    xlim([1   , 10 + 2E-2]);

ylim([0 - 3E-3 ,  1]);
set(gca,'XTick',[1:10])
% grid on
set(gca,'fontsize',fontsizes.axes);
axis normal



first_plot = true;  GOFs = NaN(length(goal_combos),length(min_GOF));
for go = 1:length(goal_combos)  %17230
    go / length(goal_combos)
    
    goals = goal_combos{go};
    Pareto_data.fullnames{goals}
    
    
    %      if any( GOFs(go,:) > [0.6 0.6 0.6] ) || all( GOFs(go,:) < [0.3 0.3 0.3] )
    %          continue
    %      end
    
    
    perfs = NaN([size(Pareto_data.interped,1) * size(Pareto_data.interped,2) , length(goals)]);
    for g  = 1:length(goals)
        perfs(:,g) = reshape( Pareto_data.interped(:,:,goals(g)) ,[],1);
    end
    
    XPc = XP(:);  YPc = YP(:);
    NaNs = any(isnan(perfs),2);
    perfs(  NaNs , :) = [];
    XPc(NaNs) = [];  YPc(NaNs) = [];
    
    %     rounded = roundn(perfs,-1);
    %     rounded = perfs;
    
    temp = round(perfs,8,'significant');
    [ p, optimal_inds] = paretoFront( temp );
    
    
    
    Optimal = NaN(size(XP));
    for oi = 1:length(optimal_inds)
        Optimal( XP == XPc(optimal_inds(oi)) & YP == YPc(optimal_inds(oi)) ) = 1;
    end
    
    A_pareto = dA * length(optimal_inds);  % approx area of Pareto region
    
    %fields = fieldnames(Pareto_data.observed); 
    fields = {'species','individuals'};
    GOF = NaN(1,length(fields));
    for ff = 1:length(fields)
        
        in_data_hull = inpolygon(XPc,YPc,Pareto_data.observed.(fields{ff}).boundary(:,1),Pareto_data.observed.(fields{ff}).boundary(:,2)); %returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
        
        in_optimal = intersect( find(in_data_hull) , optimal_inds );
        out_optimal = intersect( find(~in_data_hull) , optimal_inds);
        in_suboptimal = setdiff( find(in_data_hull) , optimal_inds );
        
        C = dA*length(in_optimal);   % area within both observations and Pareto region
        B = dA*length(out_optimal);   % area outside observations but within Pareto region
        A = dA*length(in_suboptimal); % area within observations but outside Pareto region
        
        GOF(ff) = C/(B+C) * C/(A+C);
    end
    
    GOFs(go,:) = GOF;
    
        if GOF(1) < min_GOF(1)  %all(GOF < min_GOF)
            continue
        end
    
    
    
    
    figure(52)
    if first_plot
        
        hold on
        
        pc = pcolor(XP,YP,(Optimal));
        
        
        colormap(gray(20));
        shading flat
        set(pc,'FaceAlpha',0.25);
        
%         [legend_h,object_h,plot_h,text_str] = legendflex([ individuals, species, genera, pc],  {'individuals','species','genera','Pareto-optimal'},'fontsize',12);
%         [legend_h,object_h,plot_h,text_str] = legendflex([ individuals, species, bound_individuals, bound_species, pc],  {'individual cells','species','individuals boundary','species boundary','Pareto-optimal'},'fontsize',fontsizes.legend);
         [legend_h,object_h,plot_h,text_str] = legendflex([ individuals, species,  bound_species, pc],  {'individual cells','species','species boundary','Pareto-optimal'},'fontsize',fontsizes.legend);
        
        
        set( findobj(object_h,'tag','Pareto-optimal') , 'FaceAlpha' , 0.2)
        %             set(legend_h,'Position',[153.16 677.3 266 132.5])
        hold off
        first_plot = false;
    else
        pc.CData = Optimal;
    end
    
    
    
     title({['task indices  ',num2str(goals)] , ['Goodness of Fit = ',num2str(GOF)]},'fontsize',13);
    
    drawnow
    %     if GOF > 0.5
    
    name = [];
    for tt = 1:length(goals)
        name = [name, Pareto_data.fullnames{goals(tt)}];
        if tt < length(goals)
            name = [name, ' - '];
        end
    end
%     pause
        print('-dpng','-r300',['C:\Hull\pareto\new MJ all data\','GOF ',num2str(GOF*100,'%4.0f'),'   ',num2str(length(goals)),' tasks','   ',name,'.png']);  %  curved tail
    
end