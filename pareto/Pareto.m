
Pareto_folder = 'C:\Hull\pareto\final\';

if ~exist(Pareto_folder,'dir')
    mkdir(Pareto_folder);
end

panel_goals = {[5 6 22],... eff, SNR, weighted abs gauss curv
    ... [5 6 11],... eff, SNR, unweighted abs gauss curv
    [3 6 22],... tumbling, SNR, weighted abs gauss curv (unweighted is same)
    [3 5 6 22],... eff, SNR, tumbling, weighted abs gauss curv (unweighted is same)
    [1 5 22],... uptake, tumbling, weighted abs gauss curv    (adding SNR is same, most/all constructions are same)
    [1 3 22],... uptake, eff, weighted abs gauss curv (most/all constructions are same)
    };
%%
% min_GOF = [0. 0. 0.];
min_GOF = [0.5 0.0]; %based on [ species individuals] boundaries

max_goals = 5;

fig_type = 'pareto';  % either 'survey' or 'pareto'
headliner = true;  % to generate necessary variables, run precompute_Pareto for normal resolution, then smooth_pareto.m, then precompute_Pareto.m again with high resolution
survey_SI = false; % only used for fig_type = 'survey'

panels = false;
fignum = 80;

switch fig_type
    case 'survey'
        label_selected = true; % label selected species and plot their markers in a different color
        draw_Pareto = false;
        plot_median = true;
        plot_ellipse = true;  hardcode_ellipse = true;
        do_size_scaling = true;  %plot species dots all the same size or size scaled as equivalent spherical diameter
        draw_boundary = false;
        disp_GoF = false;
        

    case 'pareto'
        label_selected = false;  % label selected species and plot their markers in a different color
        draw_Pareto = true;
        plot_median = false;
        plot_ellipse = false;
        do_size_scaling = false;  %plot species dots all the same size or size scaled as equivalent spherical diameter
        if headliner
            disp_GoF = false;
        else
            disp_GoF = true;
        end
end

plot_obs = true;
% disp_GoF = true;

selected.genus =    {'Bdellovibrio' , 'Azospirillum', 'Leptospirillum','Desulfovibrio','Desulfotomaculum', 'Vibrio'      , 'Thioalkalivibrio'  ,     'Helicobacter', 'Desulfovibrio' ,'Chlorobium'    ,  'Ammonifex' ,'Rhodospirillum'  ,     'Vibrio'   , 'Vibrio',           'Heliobacterium'   };
selected.species =  {'exovorus' ,  'halopraeferens',   'ferrooxidans', 'africanus',       'acetoxidans' ,     'ruber'          ,  'jannaschii'      ,  'flexispira'  ,    'gigas'    ,'phaeovibrioides' ,'degensii','rubrum','vulnificus' ,  'cholerae',        'undosum'  };

if strcmp(fig_type,'pareto')
if headliner
    color_Pareto = true;
    smooth_Pareto = true;  % instead of using raw Optimal matrix, completely fill in smoothed optimal region from smooth_pareto.m
    draw_boundary = false;
    do_GOF_legend = false;
else
    color_Pareto = false;
    smooth_Pareto = false;  % instead of using raw Optimal matrix, completely fill in smoothed optimal region from smooth_pareto.m
    draw_boundary = true;
    do_GOF_legend = true;
end
end
fill_straights = false;  % fill or color a thickened horizontal strip at bottom for straight rods (can be hard to see that all short straights are optimal)
if strcmp(fig_type,'survey') && survey_SI
    % for SI survey fig
    draw_boundary = true;  do_size_scaling = false; plot_ellipse = false; plot_median = false; label_selected = false; cutoff_SF1 = false;
else
    cutoff_SF1 = true;
end
% do_GOF_legend = true;

task_colors = [0 1 0; 1 0 0; 0 0 1];   %  efficiency, temporal SNR, construction ease     in_optimal;  out_optimal;   in_suboptimal;
% task_colors = [0 0 1; 1 0 0; 0 1 0];
%  task_colors = [0 1 0; 0 0 1; 1 0 0];
task_colors = [0 0 1; 0 1 0; 1 0 0]; % alternate

% task_colors = [1 0 0; 0 1 0; 0 0 1];
brighten_coeff = 0.3; %3
normalize_colors = true;
%   normalize_colors = false;
% task_colors = [repmat(0.75,2,3); 1 1 1];


if panels
    fontsizes.labels = 17;
    fontsizes.axes = 16;
    fontsizes.title = 15;
    fontsizes.caxis = 26;
    fontsizes.legend = 21;
    panel_labels = {'A','B','C','D','E'};
    
    
else
    
    fontsizes.labels = 21;
    fontsizes.axes = 21;
    fontsizes.title = 15;
    fontsizes.caxis = 26;
    fontsizes.legend = 21;
end


Pareto_h = figure(fignum);
clf; clear pc
Pareto_h.Color = 'w';



if panels
    
    precomputed_goals = {Pareto_precomputed.goals};  goal_inds = NaN(1,length(panel_goals));
    for pp = 1:length(panel_goals)
        
        for in = 1:length(Pareto_precomputed)
            if isequal( sort(panel_goals{pp}) , sort(Pareto_precomputed(in).goals))
                goal_inds(pp) = in;
                break
            end
        end
    end
end

if panels
    pans = 1:5;
else
    pans = 1;
end


for pan = pans
    
    
    if panels
        clear pc
        subtightplot(2,3,pan,[0.002 0.02],[0.06  0.01],[0.1  0.01]);
    end
    %   gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
    %            is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
    %            relatively large axis.
    %   marg_h  margins in height in normalized units (0...1)
    %            or [lower uppper] for different lower and upper margins
    %   marg_w  margins in width in normalized units (0...1)
    %            or [left right] for different left and right margins
    
    if plot_obs
        if do_size_scaling
            scale_factor = 20;  %NaN_size = NaN;
        else
            scale_factor = 1;
        end
        tmp = vertcat(data_filtered.good.(aggregation_method));     sizes = (vertcat(tmp.sph_dia)).^2;   median_size = nanmedian(sizes);
        NaN_size = median_size;
        individuals = plot([data.SF1],[data.SF2],'.','markerfacecolor',repmat(0.,1,3),'markeredgecolor',repmat(0.,1,3),'markersize', 5,'visible','off');
        hold on
        %     temp = vertcat(species_data.mean_unweighted);
        %     species = plot(temp(:,1),temp(:,2),'o','markerfacecolor','b','markersize',7);
        tmp = vertcat(data_filtered.manually_removed.(aggregation_method));
        if ~isempty(tmp)
            temp = vertcat(tmp.SF);  sizes = (vertcat(tmp.sph_dia)).^2;  sizes(isnan(sizes)) =  NaN_factor;
            %         manually_removed = plot(temp(:,1),temp(:,2),'ks','markerfacecolor','w','markersize',7);
            manually_removed = scatter(temp(:,1),temp(:,2),sizes*scale_factor,'ks','markerfacecolor','w');
        end
        %     tmp = vertcat(data_filtered.unknown_species.(aggregation_method));   temp = vertcat(tmp.SF); sizes = (vertcat(tmp.sph_dia)).^2; sizes(isnan(sizes)) =  NaN_size;
        % %     unknown_species = plot(temp(:,1),temp(:,2),'r^','markerfacecolor','w','markersize',7);
        %      unknown_species = scatter(temp(:,1),temp(:,2),sizes*scale_factor,'r^','markerfacecolor','w');
        
        tmp = vertcat(data_filtered.good.(aggregation_method));   temp = vertcat(tmp.SF);
        if do_size_scaling
            sizes = (vertcat(tmp.sph_dia)).^2;  %sizes(isnan(sizes)) =  NaN_size; % simply leave NaN size as NaN to not plot them here, instead plot separately after
        else
            sizes = repmat(35,size(temp,1),1);
        end
        %     species = plot(temp(:,1),temp(:,2),'bo','markerfacecolor','b','markersize',7);
        if panels
            scale_factor = 0.4;
        end
        species = scatter(temp(:,1),temp(:,2),sizes*scale_factor,'ko','markerfacecolor','k');
        temp2 = temp(isnan(sizes),:); sizes2 = repmat(NaN_size,sum(isnan(sizes)),1);
        species_unknown_size = scatter(temp2(:,1),temp2(:,2),sizes2*scale_factor,'ko');
        if ~isempty(data_filtered.peeled)
            temp = vertcat(data_filtered.peeled.(aggregation_method));     peeled = plot(temp(:,1),temp(:,2),'bo','markerfacecolor','w','markersize',7);
        end
        
        if label_selected
            ic = 0; clear selected_markers
            for ii = 1:length(data_filtered.good)
                if ismember(data_filtered.good(ii).genus , selected.genus) && ismember(data_filtered.good(ii).species , selected.species)
                    ic = ic + 1;
                    selected_markers(ic) = scatter(data_filtered.good(ii).(aggregation_method).SF(1), data_filtered.good(ii).(aggregation_method).SF(2), sizes(ii)*scale_factor,...
                        'marker','o','markerfacecolor',[0 0.8 0],'markeredgecolor',[0 0 0]);
                end
            end
            
            
            for ss = 1:length(selected.genus)
                for s = 1:length(data_filtered.good)
                    if strcmp( [data_filtered.good(s).genus,' ',data_filtered.good(s).species] , [selected.genus{ss},' ',selected.species{ss}])
                        selected.species_ind(ss) = s;
                        break
                    end
                end
                selected.data_inds{ss} = [];
                for s = 1:length(data)
                    if strcmp( [data(s).genus,' ',data(s).species] , [selected.genus{ss},' ',selected.species{ss}])
                        selected.data_inds{ss}(end+1) = s;
                    end
                end
                selected.fullname{ss} = [selected.genus{ss} , ' ',selected.species{ss}];
            end
            temp = vertcat(data_filtered.good(selected.species_ind).median_unweighted);   SF = vertcat(temp.SF);
            
            %         txtsel = text(SF(:,1)+0.1,SF(:,2),selected.fullname,'fontweight','bold','backgroundcolor','w','margin',eps,'fontsize',10,'fontangle','italic');
            
            
        end
        
        if plot_median
            temp = vertcat(data_filtered.good.median_unweighted);   SF = vertcat(temp.SF);  median_pt = median(SF);
            mp = scatter(median_pt(1),median_pt(2),median_size*scale_factor,'o','markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0]);
             mp2 = scatter(median_pt(1),median_pt(2),median_size*scale_factor*8,'^','markerfacecolor','none','markeredgecolor',[1 0 0]);
        end
        
        if plot_ellipse
            
            if hardcode_ellipse
%                 best_ell = [1.2839     0.097825       3.6251      0.14979       3.1327];
                best_ell = [1.3106     0.095544       3.5785      0.14945       3.1288];  % pretty confident in this one, ran overnight and met below convergence criteria
            else
                
                nvars = size(Pareto_data.observed.species.points,1);
                
                min_count = 0.5 * nvars; % region must contain at least this fraction of total # species
                % ell =      [radius_x, radius_y, center_x, center_y, angle]
                absdiffs_tol = [0.01   0.01      0.01      0.01     0.1*pi/180];
                
                nonlcon = @(ell)ell_con(ell, Pareto_data.observed.species.points, min_count);
                options = optimoptions('patternsearch','display','off');
                % guess = [0.2 0.1 5.5 0.05 45*pi/180]';
                lb = [0.5 0.05 1 0.05 0]';  ub = [5 0.25 5 0.5 2*pi]';
                n_guesses = 50;
                n_guesses = 16;
                
                best_ell = NaN(1,5);  best_fval = Inf;
                
                while true
                    
                    
                    clear ells fvals flags
%                     ppm = ParforProgressStarter2('ellipse calc', n_guesses);
                    parfor n = 1:n_guesses
                        guess = lb + (ub - lb).*rand(5,1);
                        [ells(n,:),fvals(n),flags(n)] = patternsearch(@ellobj,guess,[],[],[],[],lb,ub,nonlcon,options);
%                         if flags(n) <= 0
%                             ells(n,:) = NaN; fvals(n) = NaN;
%                         end
%                         ppm.increment(n);
                    end
%                     delete(ppm);
                    
                    ells = ells(flags > 0 , :);  fvals = fvals(flags > 0);
                    if isempty(fvals)
                        continue
                    end
                    
                    [~,inds] = sort(fvals);
                    ell = ells(inds(1),:);  fval = fvals(inds(1));
                    [  ells(inds,:)  fvals(inds)'  ];
                    [~,~,is_inside] = nonlcon(ell);
                    [min_count sum(is_inside)];
                    
%                     [fval best_fval]
                    
                    
                    if fval < best_fval
                        norm_rad = mod(best_ell(5) - ell(5) , pi);
                        angle_abs_diff = min( pi - norm_rad, norm_rad);
                        
%                         angle_abs_diff = abs( ( pi - abs(abs(best_ell(5) - ell(5)) - pi)   ) ) ;
                        
                        absdiffs = [ abs( (best_ell(1:4) - ell(1:4)) )   angle_abs_diff  ]
                        
                        
                        best_ell = ell;
                        best_fval = fval;
                    end
                    
                    
                    
                    if all( absdiffs < absdiffs_tol)
                        break
                    end
                    
                end
                
            end
            % ell = [radius_x, radius_y, center_x, center_y, angle]
            
            % https://math.stackexchange.com/questions/941490/whats-the-parametric-equation-for-the-general-form-of-an-ellipse-rotated-by-any
            t = linspace(0,1,200);
            x = best_ell(3) +    best_ell(1)* cos(2*pi*t) *cos(best_ell(5)) - best_ell(2)* sin(2*pi*t) *sin(best_ell(5));
            y = best_ell(4) +    best_ell(1)* cos(2*pi*t) *sin(best_ell(5)) + best_ell(2) *sin(2*pi*t) *cos(best_ell(5));
            hold on
            pl = plot(x,y,'r-','linewidth',0.5);
            
            
            
        end
        
        
        
        
        
        
        
        hold on
        %     temp = vertcat(genus_data.mean_unweighted);
        %     genera = plot(temp(:,1),temp(:,2),'ro','markerfacecolor','r','markersize',7);
        
        %      bound_individuals = plot(Pareto_data.observed.individuals.boundary(:,1),Pareto_data.observed.individuals.boundary(:,2),':','linewidth',1.5,'color',repmat(0.4,1,3));
        if draw_boundary
            clear bound_species
            if strcmp(fig_type,'survey')
                bound_patch = patch(Pareto_data.observed.species.boundary.cutoff_SF1(:,1),Pareto_data.observed.species.boundary.cutoff_SF1(:,2),'b');
                bound_patch.FaceAlpha = 0.2;  bound_patch.EdgeColor = 'none';
            end
            %         bound_species(1) = plot(Pareto_data.observed.species.boundary.cutoff_SF1(:,1),Pareto_data.observed.species.boundary.cutoff_SF1(:,2),'k-','linewidth',3);
            switch fig_type
                case 'survey'
                    linecolor = 'b';
                case 'pareto'
                    linecolor = 'k';
            end
            if ~panels
                linewidth =1.5;
            else
                linewidth = 0.5;
            end
            except_cutoff_line = Pareto_data.observed.species.boundary.cutoff_SF1([3:end 2],:);  cutoff_line = Pareto_data.observed.species.boundary.cutoff_SF1([2:3],:);
            bound_species(1) = plot(except_cutoff_line(:,1),except_cutoff_line(:,2),[linecolor,'-'],'linewidth',linewidth);
            bound_species(3) = plot(cutoff_line(:,1),cutoff_line(:,2),[linecolor,'--'],'linewidth',linewidth);
            bound_species(2) = plot(Pareto_data.observed.species.boundary.all_data(:,1),Pareto_data.observed.species.boundary.all_data(:,2),[linecolor,'-'],'linewidth',linewidth);
            
        end
        %     bound_genera = plot(Pareto_data.observed.genera.boundary(:,1),Pareto_data.observed.genera.boundary(:,2),'r-','linewidth',2);
        if ~panels || pan == 4
            
            xlabel('$$ \mathcal{L} $$ (elongation)','Interpreter','latex','fontsize',fontsizes.labels)
            ylabel('$$ \mathcal{K} $$ (curvature)','Interpreter','latex','fontsize',fontsizes.labels)
        end
        
        if ~panels
            set(gcf,'Position',[   467         154        1022         730]);  % aspect ratio = 1.4, but mainly need to set for axis
            set(gca,'Position',[ 0.13      0.13411        0.79      0.79]); % once fig is right aspect ratio, equal axes normalized width, height should make axes match fig aspect ratio
        else
            set(gcf,'Position',[  194    32   775   934]);
        end
        
 
        if cutoff_SF1
            set(gca,'PlotBoxAspectRatio',[   1      0.80206      0.80206]); % ALSO needed to reproduce aspect ratio...
%              pbaspect([1      0.3      0.5]);  % make aspect ratio same as performance landscape figures
       
        else
            set(gca,'PlotBoxAspectRatio',[   1      0.80206 /( max(species.XData) / 10 )    0.80206]); % ALSO needed to reproduce aspect ratio....
        end
        %         ptemp = patch('Faces',[1 2 3 4],'Vertices',[1 0.9; 1 1; 2 1; 2 0.9],'visible','off');
        hold off
        
    end
    
    switch fig_type
        case 'survey'
            if survey_SI
            grid on
            end
    end
    
    
    
    main_ax = gca;
    
    %     xlim([1   , 10 + 2E-2]);
    if cutoff_SF1
        xlim([1   , 10]);
        set(gca,'XTick',[1:10])
    else
        xlim([1 Inf]);
        set(gca,'XTick',[1 5:5:20]);
    end
    % ylim([0 - 0.01 ,  1]);
    ylim([0  ,  1]);
  
    
    % grid on
    set(gca,'fontsize',fontsizes.axes);
    if cutoff_SF1 && ~ strcmp(fig_type, 'survey')
        axis normal
    end
    % set(Pareto_h,'visible','off');
    
    if panels && pan ~= 4
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    end
    
    first_plot = true;
    if panels
        goals_loop = panel_goals(pan);
    else
        goals_loop = goal_combos;
        GOFs = NaN(length(goal_combos),length(min_GOF));
    end
    
    
    
    
    for go = 1:length(goals_loop)  %776 is the last 3-task combo
        
        
        if ~panels && length(goal_combos{go}) > max_goals
            continue
        end
        
        
        
        
        if ~panels
            go / length(goal_combos)
        end
        
        
        
        goals = goals_loop{go};
        
        if ~isequal(goals,[5 6 18])
            %                 continue
        end
        
        
        ism = ismember(construction_inds, goals);
        if sum(ism) > 1
            continue
        end
        
        if ismember(2,goals)  % fore-aft
            %         continue
        end
       
        
        %     Pareto_data.fullnames{goals};
        if ~smooth_Pareto
            perfs = NaN([size(Pareto_data.interped,1) * size(Pareto_data.interped,2) , length(goals)]);
            for g  = 1:length(goals)
                perfs(:,g) = reshape( Pareto_data.interped(:,:,goals(g)) ,[],1);
            end
        else
            perfs = NaN([size(Pareto_data_refined.interped,1) * size(Pareto_data_refined.interped,2) , length(goals)]);
            for g  = 1:length(goals)
                perfs(:,g) = reshape( Pareto_data_refined.interped(:,:,goals(g)) ,[],1);
            end
        end
        
        if smooth_Pareto
            XP = Pareto_data_refined.XP;  YP = Pareto_data_refined.YP;
        end
        
        
        XPc = XP(:);  YPc = YP(:);
        
        NaNs = any(isnan(perfs),2);
        perfs(  NaNs , :) = [];
        XPc(NaNs) = [];  YPc(NaNs) = [];
        
        %     rounded = roundn(perfs,-1);
        %     rounded = perfs;
        
        kept_inds = find(~NaNs); % indices into big original XP, YP, Optimal
        
        fields = {'species','individuals'};
        GOF = NaN(1,length(fields));
        if ~smooth_Pareto
            
            if ~panels
                ind = go;
            else
                ind = goal_inds(pan);
            end
            
            optimal_inds = Pareto_precomputed(ind).optimal_inds;
            Optimal = Pareto_precomputed(ind).Optimal;
            A_pareto = Pareto_precomputed(ind).A_pareto;  % approx area of Pareto region
            
            %fields = fieldnames(Pareto_data.observed);
            
            
            for ff = 1:length(fields)
                
                in_data_hull = inpolygon(XPc,YPc,Pareto_data.observed.(fields{ff}).boundary.cutoff_SF1(:,1),Pareto_data.observed.(fields{ff}).boundary.cutoff_SF1(:,2)); %returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
                
                in_optimal = intersect( find(in_data_hull) , optimal_inds );
                out_optimal = intersect( find(~in_data_hull) , optimal_inds);
                in_suboptimal = setdiff( find(in_data_hull) , optimal_inds );
                if ff == 1
                    region_inds = {kept_inds(in_optimal), kept_inds(out_optimal) , kept_inds(in_suboptimal)};
                end
                C = dA*length(in_optimal);   % area within both observations and Pareto region
                B = dA*length(out_optimal);   % area outside observations but within Pareto region
                A = dA*length(in_suboptimal); % area within observations but outside Pareto region
                
                GOF(ff) = C/(B+C) * C/(A+C);
            end
            
        else
            if exist('inpoly','var')  
                temp = inpoly;
            else
                temp = inpolygon(XP(:),YP(:),[main_poly.Vertices(:,1); NaN; island_poly.Vertices(:,1)],[main_poly.Vertices(:,2); NaN; island_poly.Vertices(:,2)]);
                inpoly = temp;
            end
            Optimal = NaN(size(XP));
            Optimal(temp) = 1;
        end
        
        if ~panels
            GOFs(go,:) = GOF;
        end
        
        if disp_GoF
            if first_plot
                if ~panels
                try, delete(te), end
                end
                clear text
                te = text(7.43,0.71,['GoF = ',num2str(GOF(1),2)],'fontsize',16);
                if panels
                pl = text(1.1,0.96,panel_labels{pan},'fontsize',18);
                end
                
            else
                te.String = ['GoF = ',num2str(GOF(1),3)];
            end
        end
        
        if GOF(1) < min_GOF(1)  %all(GOF < min_GOF)
            continue
        end
        
     
        set(0, 'CurrentFigure', Pareto_h);
        
        
        colors = repmat( NaN(size(Optimal)) , 1 , 1 , 3);
        
        for c = 1:3
            for i = 1:3
                temp = colors(:,:,i);
                
                temp(region_inds{c}) = task_colors(c,i);
                colors(:,:,i) = temp;
            end
        end
        
        
        if first_plot
            
            hold on
            
            
            if draw_Pareto
                
                %         pc = pcolor(XP,YP,(Optimal));
                pc = pcolor(XP,YP,Optimal);
                pc.CData = colors;
                pc.FaceAlpha = 0.5;
                if cutoff_SF1
                    pbaspect([1      0.80206      0.80206]);  % make aspect ratio same as performance landscape figures
                
                else
                    pbaspect([1      0.80206 /( max(species.XData) / 10 )     0.80206]);  % make aspect ratio same as performance landscape figures
                end
                %         colormap(gray(3));
                shading flat
                %         set(pc,'FaceAlpha',0.5);
                %         [legend_h,object_h,plot_h,text_str] = legendflex([ species],  {'species medians'},'fontsize',fontsizes.legend);
            end
            
            
            if draw_Pareto && do_GOF_legend
                for c = 1:3
                    GOF_plts(c) = patch([4 6 4],[0.1 0.1 0.2],task_colors(c,:),'linestyle','none');
                end
                
                if ~panels || pan == 1
                    
                    [GOF_leg,GOF_obj,GOF_plt,GOF_txt] = legend(GOF_plts,{'observed, optimal (R_C)','not observed, optimal (R_B)','observed, suboptimal (R_A)'},'location','northeast','fontsize',13);
                    if panels
                        set(GOF_leg, 'Position',[0.285979408969546   0.905059102008524   0.249032252603962   0.076017128515040]);
                    end
                    set( findobj(GOF_obj,'type','patch')   , 'FaceAlpha' , pc.FaceAlpha);
                end
                set(GOF_plts,'visible','off');
            end
            %         [legend_h,object_h,plot_h,text_str] = legendflex([ individuals, species, genera, pc],  {'individuals','species','genera','Pareto-optimal'},'fontsize',12);
            %         [legend_h,object_h,plot_h,text_str] = legendflex([ individuals, species, bound_individuals, bound_species, pc],  {'individual cells','species','individuals boundary','species boundary','Pareto-optimal'},'fontsize',fontsizes.legend);
            %          [legend_h,object_h,plot_h,text_str] = legendflex([ individuals, species,  bound_species, pc],  {'individual cells','species','species boundary','Pareto-optimal'},'fontsize',fontsizes.legend);
            
            %         set( findobj(object_h,'tag','Pareto-optimal') , 'FaceAlpha' , 0.2)
            %             set(legend_h,'Position',[153.16 677.3 266 132.5])
            hold off
            if ~panels
            first_plot = false;
            end
        else
            pc.CData = colors;
            
            
            %         delete(pc);
            %         hold on
            %         pc = pcolor(XP,YP,Optimal);
            %
            %         pc.CData = Optimal;
            %
            %         pc.CData = colors;
            %         pc.FaceAlpha = 0.5;
            
            %         colormap(gray(3));
            %         shading flat
            hold off
            %          [GOF_leg,GOF_obj,GOF_plt,GOF_txt] = legend(GOF_plts,{'observed, optimal (C)','not observed, optimal (B)','observed, suboptimal (A)'},'location','northeast');
        end
        
        
        
        if draw_Pareto && color_Pareto
            
            
            Optimal2 = false(size(Optimal));
            
            Optimal2(Optimal == 1) = true;
            
            Optimal_weights = NaN([size(Optimal) 3]);  Optimal_color = Optimal_weights;
            
            clear opt_metric_normalized
            
            clear metric_ranges
            for task = 1:3
                if ~smooth_Pareto
                    metric = squeeze(Pareto_data.interped(:,:,goals(task)));
                else
                    metric = squeeze(Pareto_data_refined.interped(:,:,goals(task)));
                end
                
                
                opt_metric = metric(Optimal2);
                
                metric_range = [min(opt_metric(:)) max(opt_metric(:))];
                metric_ranges(task,:) = metric_range;
                
                opt_metric_normalized = ( opt_metric - metric_range(1) ) / (metric_range(2) - metric_range(1));
                full_normalized = NaN(size(Optimal));
                full_normalized(Optimal2) = opt_metric_normalized;
                
                Optimal_weights(:,:,task) = full_normalized;
                
                
                
            end
            
            
            
            
            for nr = 1:size(Optimal_color,1)
                for nc = 1:size(Optimal_color,2)
                    Optimal_color(nr,nc,:) = sum( repmat(reshape( Optimal_weights(nr,nc,:) , [3 1 1] ),1,3) .* task_colors , 1);
                end
            end
            
            
            
            
            %%
            %         normalization_constants = [0.2989  0.5870  0.1140];
            normalization_constants = [1 1 1];
            I = 1;  % total "intensity" of each normalized color
            
            %         [r g b]
            %         c1*r + c2*g + c3*b = I
            %         m*(c1*r + c2*g + c3*b) = I
            %         normalized = [m*r m*g m*b]
            
            norm_colors = @(colors) repmat( I ./ sum( repmat(  reshape(normalization_constants,[1 1 3])  ,  size(colors,1) , size(colors,2) , 1  )  .*  colors  , 3)  , 1,1,3)  .* colors;
            
            if normalize_colors
                normalized_color = norm_colors(Optimal_color);
            else
                normalized_color = Optimal_color;
            end
            
            normalized_color = brighten(normalized_color, brighten_coeff);
            
            
            %         figure(234)
            %         for iii = 1:3
            %             subplot(1,4,iii)
            %             temp = normalized_color(:,:,iii);
            %             histogram(temp(:));
            %         end
            %         subplot(1,4,4)
            %         temp = sum( repmat(  reshape(normalization_constants,[1 1 3])  ,  size(normalized_color,1) , size(normalized_color,2) , 1  )  .*  normalized_color  , 3);
            %         histogram(temp(:));
            
            
            
            
            pc.CData = normalized_color;
            %%
            %         pc.CData = Optimal_color;
            
            
            pc.FaceAlpha = 1;
            shading interp
            
            
            
% add slightly thick line around boundary of main region to fix problem of
% region near sphere not appearing to be optimal (due to narrowness)
% https://undocumentedmatlab.com/blog/plot-line-transparency-and-color-gradient
good = ~isnan( pc.CData(:,:,1) );

clear Fc linecolor
for i = 1:3
ctemp = pc.CData(:,:,i);

Fc{i} = scatteredInterpolant(XP(good),YP(good),ctemp(good),'linear','linear');

linecolor(:,i) = Fc{i}(upper_boundary_splined(:,1),upper_boundary_splined(:,2));

end

linecolor(linecolor < 0) = 0;  linecolor(linecolor > 1) = 1;
% linecolor0 = linecolor;

linecolor = linecolor * 255;

linecolor = uint8(linecolor);

ColorData = [linecolor' ; uint8(repmat(255,1,size(linecolor,1))) ];
hold on
smoothed_main = plot(upper_boundary_splined(:,1),upper_boundary_splined(:,2),'k-','LineWidth',2);
drawnow
set(smoothed_main.Edge,'ColorBinding','interpolated','ColorData',ColorData);
   set(gca,'ClippingStyle', 'rectangle');
   drawnow
            
            
            
            %%
            inset = axes(Pareto_h,'position',[  0.65          0.7          0.15          0.15],'visible','off');
            
            %%
            Vertices = [ 0.5 sqrt(3)/2;  1 0;  0 0;];   Vertices(:,end+1) = 0; % subdivide_tri needs 3D coords
            
            Faces = [1:size(Vertices,1)];
            
            TR = triangulation(Faces,Vertices);  % needed to get barycentric coords for relative perf weights across triangle
            
            
            n = 9;  % # subdivisions of triangular mesh
            [ Vertices, Faces ] = subdivide_tri( Vertices, Faces );
            for i = 2:n
                [ Vertices, Faces ] = subdivide_tri(  Vertices, Faces  );
            end
            
            %         col=[1 0 0; 0 1 0; 0 0 1]; % RGB colors at corners for each task
            %         col = [1 0 1; 1 1 0; 0 1 1];
            col = task_colors;
            
            perf_weights = cartesianToBarycentric(TR,ones(size(Vertices,1),1),Vertices);  perf_weights(perf_weights < 0) = 0;  % fix numerical roundoff
            %         perf_weights = 0
            perf_colors = NaN(size(perf_weights));
            for i = 1:size(perf_weights,1)
                perf_colors(i,:) = sum( repmat( perf_weights(i,:)' , 1, 3) .* col , 1 );
            end
            
            I = 1;  normalization_constants = [1 1 1];
            %          normalization_constants = [0.2989  0.5870  0.1140] / 0.1140;
            %         normalized_weights = perf_weights *1;
            
            %   perf_colors = brighten(perf_colors , 0.3);
            
            %         normalized_colors = repmat( I ./ sum( repmat(  normalization_constants  ,  size(perf_colors,1) , 1   )  .*  perf_colors  , 2)  , 1, 3)  .* perf_colors;
            %         normalized_weights = normalized_weights / max(normalized_weights(:));
            
            normalized_colors =    perf_colors;
            
            normalized_colors = brighten(normalized_colors , brighten_coeff);
            
            
            %         figure(235)
            %         for iii = 1:3
            %             subplot(1,4,iii)
            %             temp = normalized_colors(:,iii);
            %             histogram(temp(:));
            %         end
            %         subplot(1,4,4)
            %         temp = sum( repmat(  normalization_constants  ,  size(normalized_colors,1) , 1   )  .*  normalized_colors  , 2);
            %         histogram(temp(:));
            
            
            % normalized_colors = normalized_colors *1.6;
            
            %          figure(723);    kp = patch( 'Vertices',Vertices, 'Faces',Faces, 'EdgeColor','none','FaceVertexCData', normalized_colors,'FaceColor','interp');  axis equal; axis tight
            %%
            figure(Pareto_h);
            kp = patch(inset, 'Vertices',Vertices, 'Faces',Faces, 'EdgeColor','none','FaceVertexCData', normalized_colors,'FaceColor','interp');
            axis equal; axis tight
            clear text
            xy=[0 0; 1 0; 0.5 sqrt(3)/2];
            triangle_labels(1) = text(xy(1,1) - 0.1, xy(1,2)-0.2 + 0.0,{'Construction', 'Ease'},'HorizontalAlignment', 'center','FontSize',fontsizes.legend - 7);
            triangle_labels(2) = text(xy(2,1) + 0.1, xy(2,2)-0.2 + 0.0,{'Chemotactic', 'SNR'},'HorizontalAlignment', 'center','FontSize',fontsizes.legend - 7);
            triangle_labels(3) = text(xy(3,1), xy(3,2)+0.1 + 0.05,{'Swimming Efficiency'},'HorizontalAlignment', 'center','FontSize',fontsizes.legend - 7);
            kp.FaceAlpha = 1;
            
            
            %         rect = annotation('rectangle',[0.67786 (0.45575 + 0.02) 0.22 0.36318],'LineWidth',2);
            
            %         rect_title = annotation('textbox','linestyle','none','fontweight','bold','fontsize',16,...
            %             'String',{'Relative weight per task'},'Position',[0.70853 (0.72551 + 0.02) 0.3 0.085561]);
            
            title('Performance tradeoffs between Pareto-optimal morphotypes','fontsize',fontsizes.title);
        end
        
        if fill_straights
            short_straights = pc.XData <= Inf & pc.YData == 0;
            x = pc.XData(short_straights);
            y = pc.YData(short_straights);
            z = zeros(size(x));
            
            if draw_Pareto && color_Pareto
                clear colors
                for i = 1:3
                    ctemp = squeeze(pc.CData(:,:,i));
                    color = ctemp(short_straights);
                    colors(:,:,i) = [color'; color'];
                end
                
                
                
                short_straight_segment = surface(main_ax,[x'; x'],[y'; y'],[z'; z'],colors,...
                    'facecolor','flat',...
                    'edgecolor','interp',...
                    'linew',5,'facealpha',1,'edgealpha',1);
            elseif draw_Pareto
                short_straight_segment = surface(main_ax,[x'; x'],[y'; y'],[z'; z'],ones(size([x'; x'])),...
                    'linew',5,'facealpha',1,'edgealpha',1);
                color = repmat(0.75,1,3);
                short_straight_segment.CData = repmat( reshape(color,[1 1 3]) , [size([x'; x']) 1]);
                short_straight_segment.EdgeColor = color;
                
            end
        end
        
        
        uistack(species,'top');
        if draw_boundary
            uistack(bound_species,'top');
        end
        if label_selected
            uistack(selected_markers,'top'); % make sure we see selected markers and not original markers
        end
        if ~draw_Pareto
            pc.Visible = 'off';
        end
        
        %      title({['task indices  ',num2str(goals)] , ['Goodness of Fit = ',num2str(GOF)]},'fontsize',13);
        
        drawnow
        %     if GOF > 0.5
%         axes(main_ax);
        
        name = [];
        for tt = 1:length(goals)
            name = [name, Pareto_data.fullnames{goals(tt)}];
            if tt < length(goals)
                name = [name, ' - '];
            end
        end
        name
%                               pause
        if strcmp(fig_type,'pareto')
%             print('-dpng','-r200',[Pareto_folder,'GOF ',num2str(GOF(1)*100,'%4.0f'),'   ',num2str(length(goals)),' tasks','   ',name,'.png']);  %  curved tail
        end
        if strcmp(fig_type,  'survey')
            break
        end
    end
    
    
end
% this seems to work well, though it's not vector
% export_fig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\Pareto headliner.png','-r600','-nocrop')

%  export_fig(gcf,'C:\allsync\all papers\curved rod ms\figures\temp\GoF comparisons.png','-r600','-nocrop')
 
 