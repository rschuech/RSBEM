


%%
nAR1 = 150;  nAR2 = 150;
 nAR1 = 150;  nAR2 = 150;
dA = (10 - 1)/(nAR1 - 1) * (1 - 0)/(nAR2 - 1);  %area of a rectangle centered on each interrogation point

[XP,YP] = ndgrid(linspace(1,10,nAR1),linspace(0,1,nAR2));
[standardized_XPYP] = standardize([XP(:) YP(:)], [1 10; 0 1]);
standardized_XP = NaN(size(XP));  standardized_XP(:) = standardized_XPYP(:,1);
standardized_YP = NaN(size(YP));  standardized_YP(:) = standardized_XPYP(:,2);


clear interped
fields = fieldnames(Pareto_data);
for f = 1:length(fields)
    field = fields{f};
    %     [~,I] = sort( Pareto_data.(field).Z(:) ,'ascend' );
    interped.(field) = Pareto_data.(field).F(standardized_XP,standardized_YP);
    [~,I] = sort( interped.(field)(:) ,'ascend' );
    %     [~,I] = sort( Pareto_data.(field).depvar ,'ascend' );
    inds = 1:length(I);
    %     Pareto_data.(field).rank = NaN(size( Pareto_data.(field).Z ));
    %     Pareto_data.(field).rank = NaN(size( Pareto_data.(field).depvar ));
    Pareto_data.(field).rank = NaN(size( interped.(field) ));
    Pareto_data.(field).rank(I) = inds;  %I rearranges inds so that we get the sorted ranks
    %     Pareto_data.(field).rank(isnan(Pareto_data.(field).Z)) = NaN;
    Pareto_data.(field).rank(isnan(interped.(field))) = NaN;
end
%%



goals = {'Construction_Ease','Power_eff','Dm','temporal_SN','fore_aft_SN', 'tumbling'  ,  'uptake'};
goals_titles = {'Construction Ease','Swimming Efficiency','Dispersal','Temporal S/N','Fore-Aft S/N','Tumbling Ease'  ,  'Uptake'};
names_titles = {'Construction','Swimming','Dispersal','Temporal','Fore-Aft','Tumbling' , 'Uptake'};


goals = {'Construction_Ease',  'uptake', 'Dbr'};
goals_titles = {'Construction Ease' ,  'Uptake','Brownian Dispersal'};
names_titles = {'Construction', 'Uptake','Dispersal'};

goals = {'Construction_Ease',  'test1', 'test2'};
goals_titles = {'Construction Ease' ,  'test1','test2'};
names_titles = {'Construction', 'test1','test2'};

goals = {'Construction_Ease',  'test1', 'uptake'};
goals_titles = {'Construction Ease' ,  'test1','uptake'};
names_titles = {'Construction', 'test1','uptake'};

goals = {'Construction_Ease',  'Dbr', 'test2'};
goals_titles = {'Construction Ease' ,  'Dbr','test2'};
names_titles = {'Construction', 'Dbr','test2'};

% goals = {'translation_ease','rotation_resistance'};
% goals_titles = {'translation_ease','rotation_resistance'};
% names_tiles = {'translation_ease','rotation_resistance'};

% goals = {'Construction_Ease','Power_eff','fore_aft_SN'};
% goals_titles = {'Construction Ease','Swimming Efficiency','Fore-Aft S/N'};
% names_titles = {'Construction','Swimming','Fore-Aft'};

% goals = {'Construction_Ease1','Construction_Ease2','Power_eff','Dm','temporal_SN','fore_aft_SN', 'tumbling'  ,  'uptake'};
% goals_titles = {'Construction Ease1','Construction Ease2','Swimming Efficiency','Dispersal','Temporal S/N','Fore-Aft S/N','Tumbling Ease'  ,  'Uptake'};
% names_titles = {'Construction1','Construction2','Swimming','Dispersal','Temporal','Fore-Aft','Tumbling' , 'Uptake'};

goal_combos = {};  title_combos = {};  name_combos = {};
for ngoals = length(goals):-1:3
    temp = combnk(goals,ngoals);
    temp2 = combnk(goals_titles,ngoals);
    temp3 = combnk(names_titles,ngoals);
    
    for t = 1:size(temp,1)
        
        goal_combos = [goal_combos; {temp(t,:)}];
        tit = []; name = [];
        for tt = 1:length(temp2(t,:))
            tit = [tit , temp2{t,tt}];
            name = [name,temp3{t,tt}];
            if tt < length(temp2(t,:))
                tit = [tit, '  ,  '];
                name = [name, '_'];
            end
        end
        title_combos = [title_combos; {tit}];
        name_combos = [name_combos; {name}];
    end
end


% goals = {'Construction_Ease','Power_eff','Dm'};
% goals = {'Construction_Ease','Power_eff','Dm','temporal_SN','fore_aft_SN'};
%goals = {'temporal_SN','Power_eff','fore_aft_SN','Dm','tumbling'};

for go = 1:length(goal_combos)
    go / length(goal_combos)
    goals = goal_combos{go};
    
    
    perfs = [];
    for g  = 1:length(goals)
        perfs(:,g) = interped.(goals{g})(:);
    end
    
    XPc = XP(:);  YPc = YP(:);
    NaNs = any(isnan(perfs),2);
    perfs(  NaNs , :) = [];
    XPc(NaNs) = [];  YPc(NaNs) = [];
    
    %     rounded = roundn(perfs,-1);
    rounded = perfs;
    
    
    [ p, optimal_inds] = paretoFront(rounded );
    
    % remove points with buckling over the assumed hard limit from the optimal
    % set
    
    %  optimal_inds = setdiff(optimal_inds,buckling_too_high);
    
    
    Optimal = NaN(size(XP));
    for oi = 1:length(optimal_inds)
        Optimal( XP == XPc(optimal_inds(oi)) & YP == YPc(optimal_inds(oi)) ) = 1;
    end
    
    A_pareto = dA * length(optimal_inds);  % approx area of Pareto region
    
    in_data_hull = inpolygon(XPc,YPc,observed.individuals(:,1),observed.individuals(:,2)); %returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    in_data_hull = inpolygon(XPc,YPc,observed.species(:,1),observed.species(:,2)); %returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    
    in_optimal = intersect( find(in_data_hull) , optimal_inds );
    out_optimal = intersect( find(~in_data_hull) , optimal_inds);
    in_suboptimal = setdiff( find(in_data_hull) , optimal_inds );
    
    C = dA*length(in_optimal);
    B = dA*length(out_optimal);
    A = dA*length(in_suboptimal);
    
    GOF = C/(B+C) * C/(A+C);
    
    GOFs(go) = GOF;
    
    if GOF < 0.0
        continue
    end
    
    %%
    
    measurements_only = false;
    plot_obs = false;
   
    
    fontsizes.labels = 18;
    fontsizes.axes = 18;
    fontsizes.title = 15;
    fontsizes.caxis = 26;
    fontsizes.legend = 22;
    
    if ~measurements_only
        figure(825)
    else
        figure(830)
    end
    clf
    
    
    if ~measurements_only
        % [Contour_data,contour_handle] = contour(unique(X),unique(Y),Pareto_data.fore_aft_SN.Z,[1 1],'k--','linewidth',1);
        Contour_data = contourcs(unique(X),unique(Y),Pareto_data.fore_aft_SN.Z,[1 1]);
        % contour_handle = plot(Contour_data(2).X,Contour_data(2).Y,'k--','linewidth',1);
        
        
        hold on
        
        
        pc = pcolor(XP,YP,(Optimal));
        
        set(pc,'FaceAlpha',0.15);
        colormap(gray(20));
        shading flat
        
    end
    if plot_obs
    %        po = plot(XPc(out_optimal),YPc(out_optimal),'bo','markerfacecolor','b','markersize',8);
    temp = vertcat(plot_data.mean_unweighted);
    all_meas = plot(meas(:,1),meas(:,2),'ko','markerfacecolor','k','markersize',1);  % individual dots
    % bound_individuals = plot(observed.individuals(:,1),observed.individuals(:,2),'k--','linewidth',1.5); %individuals hull
    hold on
    species_meas = plot(temp(:,1),temp(:,2),'ko','markerfacecolor','k','markersize',10); %species circles
    if ~measurements_only
        %bound_species = plot(observed.species(:,1),observed.species(:,2),'k-','linewidth',1.5); %species hull
        bound_species = plot(observed.species(1:2,1),observed.species(1:2,2),'k-','linewidth',1.5); %species hull
        plot(observed.species(3:end,1),observed.species(3:end,2),'k-','linewidth',1.5)
    end
    
    end
    axis normal
    
    
    if ~measurements_only && plot_obs
        % plot edges of parameter space
        try, delete(hullh), end;
        xtemp = X(~isnan(Pareto_data.Power_eff.Z));  ytemp = Y(~isnan(Pareto_data.Power_eff.Z));
        k = convhull(xtemp,ytemp,'simplify',true);
        hullpts = [xtemp(k) ytemp(k)];
        hullpts(hullpts(:,1) == 10 & hullpts(:,2) <= 0.95 - 1E-2 , :) = []; % remove right vertical boundary
        hullpts(hullpts(:,2) == 0 & hullpts(:,1) > 1 , :) = [];  %remove bottom horizontal boundary
        hullpts(1,:) = [];  %get rid of annoying connecting line
        hullh = plot(hullpts(:,1),hullpts(:,2),':','linewidth',2,'color',[repmat(0.7,1,3)]);
        
        
        % infeasible region based on fore-aft S/N
%         [~,ind] = near(Contour_data(2).X(1),hullpts(:,1));
%         verts = [flipud(hullpts(1:ind,:)); [Contour_data(3).X' Contour_data(3).Y']; flipud([Contour_data(2).X' Contour_data(2).Y'])];
%         faces = 1:length(verts);
%         infeasible = patch('Faces',faces,'Vertices',verts,'visible','off');
%         %     infeasible_h = hatchfill2(infeasible,'HatchSpacing',10,'HatchLineWidth',0.5);
%         infeasible_h = hatchfill2(infeasible,'cross');
%         set(infeasible_h,'visible','on','color',repmat(0.8,1,3));
    end
    
    if ~measurements_only
        xlim([1 10+1.5E-2]);
    else
        xlim([1 26]);
    end
    ylim([0 1]);
    grid on
    set(gca,'fontsize',fontsizes.axes);
    
    %     title({title_combos{go} , ['Goodness of Fit = ',num2str(GOF)]},'fontsize',13);
    
    
    
    xlabel('AR_1 (elongation)','fontsize',fontsizes.labels)
    ylabel('AR_2 (curvature)','fontsize',fontsizes.labels)
    set(gcf,'Position',[ 680          82        1051         896]);
    if ~measurements_only
        set(gca,'XTick',[1:10])
        
        
        %     leg = legend([all_meas, species_meas, bound, buckling_handle, hullh],'individual cells','species averages','Pareto optimal region','max buckling force 25% higher than sphere','geometrically possible');
        %       leg = legend([all_meas, species_meas, bound, hullh],'individual cells','species averages','Pareto optimal region','geometrically possible');
        %     set(leg,'Position',[ 0.38503      0.71028      0.44244      0.16295]);
        ptemp = patch('Faces',[1 2 3 4],'Vertices',[1 0.9; 1 1; 2 1; 2 0.9],'visible','off');
        % [leg,icons] = legend([ all_meas, species_meas, bound_species, pc, ptemp, hullh],  {'individuals','species means','species means outline','Pareto-optimal','worse Fore-Aft S/N than sphere','geometrically possible'},'fontsize',12);
        % tex = findobj(icons,'type','text');  set(tex,'fontsize',12);
        % pat = findobj(icons,'type','patch');  set(pat,'FaceAlpha',0.15);  % 8th entry should be optimal patch
        
        [legend_h,object_h,plot_h,text_str] = legendflex([ all_meas, species_meas, bound_species, pc, ptemp, hullh],  {'individuals','species means','species means outline','Pareto-optimal','worse Fore-Aft S/N than sphere','geometrically possible'},'fontsize',12);
        set( findobj(object_h,'tag','Pareto-optimal') , 'FaceAlpha' , 0.2)
        handle = findobj(object_h,'tag','worse Fore-Aft S/N than sphere');
        set(handle,'visible','off')
        infeasible_hleg = hatchfill2(handle,'cross');
        set(infeasible_hleg,'visible','on','color',repmat(0.8,1,3));
    else
        [legend_h,object_h,plot_h,text_str] = legendflex([ all_meas, species_meas],  {'individuals','species means'},'fontsize',12);
    end
    set(legend_h,'Position',[153.16 677.3 266 132.5])
    %     hold off
    %    pause
    drawnow
        if GOF > 0.5
    print('-dpng','-r300',['E:\Hull\pareto\nonmotile\','GOF ',num2str(GOF*100,'%4.0f'),'%   ',name_combos{go},'.png']);  %  curved tail
        end
    %pause
    % shateroonie
    %  print('-dpng','-r300',['E:\Hull\graphs2\',name_combos{go},' straight opt tail','.png']);  %  curved tail
    
    % print('-dpng','-r300',['E:\Hull\pareto\','GOF ',num2str(GOF*100,'%4.0f'),'%   ',name_combos{go},'.png']);  %  curved tail
    %%
end