extrap_method = 'none';

% interp_method = 'linear';
interp_method = 'natural';

do_save = false;

folder = 'C:\Hull\Results\fixed\';

% files = {'Uptake_Interpolant.mat',...
%     'Results_Fore_Aft_SN_synced2.mat',...
%     'Results_Tumbling_Ease.mat',...
%     'Results_median_tail_Dm.mat',....mat',...          Results_Dm_coarse_redo_1_incl_coarse_redo_fixed.mat
%     'Results_median_tail_eff.mat',....mat',...       'Results_efficiency_uber_refined_sweep1_synced.mat'
%     'Results_median_tail_temporal_synced.mat',...      Results_Temporal_SN_radial_guesses_1_synced2.mat
%     'Construction_Ease_new.mat'};

files = {'Uptake_Interpolant.mat',...
    'Results_Fore_Aft_SN_synced2.mat',...
    'Results_Tumbling_Ease.mat',...
    'Results_Dm_coarse_redo_1_incl_coarse_redo_fixed.mat',....mat',...          Results_Dm_coarse_redo_1_incl_coarse_redo_fixed.mat
    'Results_efficiency_uber_refined_sweep1_synced.mat',....mat',...       'Results_efficiency_uber_refined_sweep1_synced.mat'
    'Results_Temporal_SN_radial_guesses_1_synced2.mat',...      Results_Temporal_SN_radial_guesses_1_synced2.mat
    'Construction_Ease_new.mat'};


Depvars = {'Uptake','Fore-Aft SN','Tumbling Ease','Dispersal','Swimming Efficiency', 'Temporal SN','Construction Ease'}; %,'tau_a','lambda_a','error_angle'};
Fnames = {'uptake',  'fore_aft',    'tumbling',           'Dm',                'eff',          'temporal',    'construction'  };
interp_methods = repmat({interp_method},1,length(Fnames));



clear F Best

for d = 1:length(Depvars)
    depvar = Depvars{d};
    clear Results
    load([folder,files{d}]);
    interp_method = interp_methods{d};
    
    if exist('Results','var') % is this a metric that required simulations
        
        
        for i = 1:length(Results)
            Results(i).AR1 = roundn(Results(i).AR1,-10);
            Results(i).AR2 = roundn(Results(i).AR2,-10);
        end
        
        AR1_AR2 = [ [Results.AR1]' [Results.AR2]' ];
        AR1_AR2_unq = unique(AR1_AR2,'rows');
        %         if d == 4
        %             bads = [5.5 0.55; 6.5 0.3; 7 0.25;];
        %         AR1_AR2_unq = setdiff(AR1_AR2_unq , bads,'rows');
        %         end
        
        AR1_AR2_data.(Fnames{d}) = AR1_AR2_unq;
        
        clear best
        for i = 1:size(AR1_AR2_unq,1) %each unique body shape
            body = AR1_AR2_unq(i,:);
            inds = find(AR1_AR2(:,1) == body(1) & AR1_AR2(:,2) == body(2)); %all the runs with this body shape
            Results_subset = Results(inds);
            
            switch depvar
                case 'Tumbling Ease'
                    metric = 1./[Results_subset.tau_a];
                    
                case 'Motile Diffusivity'
                    metric = [Results_subset.Dm];
                case 'Swimming Efficiency'
                    metric = [Results_subset.Power_eff];
                    %                      metric = [Results_subset.Adj_Speed];
                case 'Fore-Aft SN'
                    taxis = [Results_subset.taxis];
                    fore_aft = [taxis.fore_aft];
                    fore_aft_SN = [fore_aft.SN];
                    metric = fore_aft_SN;
                case 'Temporal SN'
                    taxis = [Results_subset.taxis];
                    temporal = [taxis.temporal];
                    temporal_SN = [temporal.SN];
                    metric = temporal_SN;
            end
            
            
            [~,ind] = max(metric);
            best.body(i,:) = body;
            best.metric(i,1) = metric(ind);
            vars = {'amp','lambda','nlambda','Avg_Omega','Adj_Speed','tau_a','avg_swimming_axis'};
            for v = 1:length(vars)
                if ~isfield(Results_subset,vars{v}) || isempty([Results_subset.(vars{v})])
                    best.(vars{v})(i,:) = NaN;
                else
                    var_subset = [Results_subset.(vars{v})]';
                    best.(vars{v})(i,:) = var_subset(ind,:);
                end
            end
            
            
        end
        Best.(Fnames{d}) = best;
        
        limits = [1 10; 0 1];
        [standardized_bodies] = standardize(best.body, limits);
        
        %         temp = round(best.metric,3,'significant');
        
        F.(Fnames{d}).F.metric = scatteredInterpolant(standardized_bodies(:,1), standardized_bodies(:,2), best.metric,interp_method,extrap_method);
        vars2 = {'amp','lambda','nlambda','Avg_Omega','Adj_Speed','tau_a'};
        for v = 1:length(vars2)
            F.(Fnames{d}).F.(vars2{v}) = scatteredInterpolant(standardized_bodies(:,1), standardized_bodies(:,2),  best.(vars2{v}),interp_method,extrap_method);
            
        end
    elseif strcmp(depvar,'Uptake')
        Uptake_Interpolant.Method = interp_method;
        %         Uptake_Interpolant.Values = round(Uptake_Interpolant.Values , 2, 'significant');
        F.(Fnames{d}).F.metric = Uptake_Interpolant;
        %         F.(Fnames{d}).F.metric = F_fix;
    elseif strcmp(depvar,'Construction Ease')
        F.(Fnames{d}) = Construction_Eases; %contains many sub-fields for different construction ease definitions
    end
    
end


%% Plot interpolants for each metric
% 0.7 0.8 0.9 0.95 0.97 0.99      1.01 1.015  1.016 1.0175 1.018 1.02
cvals = { ...
    [1.02 1.1:0.1:1.5],... uptake
    [0:1:5],... fore-aft
    [0.1:0.1:1],... tumbling
    [1.03 1.9  1.1:0.1:3],... Dm  1.1:0.1:2
    [0.7 0.8 0.9 0.95 0.97 0.99      1.01 1.015  1.016 1.0175 1.02],... efficiency  value of 1 done as special case below to highlight it
    [ 1.1:0.1:3 ],... temporal
    [],[],[],[],[],[],[],[],[],[],[],[  0.9 0.8:-0.1:0.1 ],[] , []}; % construction
%   [],[],[],[],[],[1E-4 1E-3 1E-2],[],[],[],[],[],[1E-4 1E-3 1E-2 ],[]}; % construction
%
%  [0.85 0.9 0.95 0.97 0.98 0.99 1  1.006 1.0075 1.008 1.009 ],... adj speed
fancy_names = {'Uptake','Fore-Aft SN','Tumbling Ease','Dispersal','Swimming Efficiency', 'Temporal SN',...
    {'Construction Ease','via surface-averaged mean curvature'},{'Construction Ease','via globally max curvature'},{'Construction Ease','via surface-averaged mean absolute curvature'},...
    {'Construction Ease','via surface-averaged max curvature'},{'Construction Ease','via circumferential principle curvature'},{'Construction Ease','via axial principle curvature'},...
    {'Construction Ease','via surface-averaged mean curvature','area weighted'},{'Construction Ease','via globally max curvature','area-weighted'},...
    {'Construction Ease','via surface-averaged mean absolute curvature','area-weighted'},{'Construction Ease','via surface-averaged max curvature','area-weighted'},...
    {'Construction Ease','via circumferential principle curvature','area-weighted'},{'Construction Ease','via axial principle curvature','area-weighted'},...
    {'Construction Ease','based on only area'}  ,{'Construction Ease','test'} };

fancy_names{18} = 'Construction Cost of Cell Shape';  fancy_names{6} = 'Finding Nutrients';  % overwrite title for poster

npts = 150;

fullnames = {};  metrics = {}; displaynames = {};
try, delete(handles); end
Fs = {};  handles = [];
% options = optimset('MaxFunEvals',1E9,'MaxIter',1E9,'TolX',1E-6,'TolFun',1E-8,'display','iter');
options = optimoptions('patternsearch','display','off','FunctionTolerance',1E-8,'MaxFunctionEvaluations',1E9,'MaxIterations',1E9,'MeshTolerance',1E-8,'StepTolerance',1E-8);
normalize_metrics = true; % normalize to metric value for a sphere

fields = fieldnames(F);


for f = 1:length(fields)
    metric = fields{f};
    displayname = Depvars{f};
    
    
    temp = fieldnames(F.(metric));  %isfield doesn't work for some reason....
    if ~any(ismember(temp,'F')) % are there different versions of this metric (e.g. Construction Ease)?
        fields2 = fieldnames(F.(metric));
        
        for ff = 1:length(fields2)
            metric2 = fields2{ff};
            temp2 = fieldnames(F.(metric).(metric2));
            if ~any(ismember(temp2,'F')) % are there different versions of this metric (e.g. Construction Ease)?
                fields3 = fieldnames(F.(metric).(metric2));
                
                for fff = 1:length(fields3)
                    metric3 = fields3{fff};
                    temp3 = fieldnames(F.(metric).(metric2).(metric3));
                    if ~any(ismember(temp3,'F'))
                        stopafra % will need to modify code for more nested fields
                    end
                    
                    fullnames{end+1} = [metric,'  ',metric2,'  ',metric3];
                    Fs{end+1} = F.(metric).(metric2).(metric3).F.metric;
                    metrics{end+1} = metric;
                    displaynames{end+1} = displayname;
                end
                
            else
                fullnames{end+1} = [metric,'  ',metric2];
                Fs{end+1} = F.(metric).(metric2).F.metric;
                metrics{end+1} = metric;
                displaynames{end+1} = displayname;
            end
        end
    else
        fullnames{end+1} = [metric];
        Fs{end+1} = F.(metric).F.metric;
        metrics{end+1} = metric;
        displaynames{end+1} = displayname;
        
    end
    
end




[X,Y] = meshgrid(linspace(limits(1,1),limits(1,2),npts),linspace(limits(2,1),limits(2,2),npts));
[standardized_XY] = standardize([X(:) Y(:)], limits);
standardized_X = NaN(size(X));  standardized_X(:) = standardized_XY(:,1);
standardized_Y = NaN(size(Y));  standardized_Y(:) = standardized_XY(:,2);

fontsizes.labels = 16;
clear clab axs;  axi = 0; % counter for selected optimal shape axes
fignum = 99;

for m =  [1]   % 1:length(Fs) % 5 6  18
    %     switch metrics{m}
    %         case 'uptake'
    %             limits = [1 10; 0 0.985];
    %         otherwise
    limits = [1 10; 0 1];
    %     end
    
    if ~normalize_metrics
        
        Z = Fs{m}(standardized_X,standardized_Y);
    else
        if strcmp(metrics{m},'construction')
            %             Z = log10( Fs{m}(standardized_X,standardized_Y) )/ log10( Fs{m}(0,0) );
            Z = Fs{m}(standardized_X,standardized_Y) / Fs{m}(0,0);
        else
            Z = Fs{m}(standardized_X,standardized_Y) / Fs{m}(0,0);
            %               Z = Fs{m}(standardized_X,standardized_Y) ;
        end
    end
    
    
    
    
    
    fignum = fignum + 1;
    handles(m) = figure(fignum);  clf;  set(gcf,'Position',[ 683   555   596   421]);
    set(handles(m),'Name',[ fullnames{m} ]);
    
    
    pcolor(X,Y,Z);  shading interp;  cb = colorbar; 
    if ismember(m,[6 4 1])
        cb.Label.String = 'performance relative to spherical';
    cb.Label.FontSize = fontsizes.labels + 4;
    end
    
    title(fancy_names{m},'interpreter','none','fontsize',20);
    
    
    switch m
        case {7:10 13:16}
            
            caxis([0.99 1]);
    end
    
    grid on
    hold on
    
    if ~isempty(cvals{m})
        [C,ch] = contour(unique(X),unique(Y),Z,cvals{m},'k--','linewidth',1);
        if strcmp(metrics{m},'eff')
        
            [C_temp,ch_temp] = contour(unique(X),unique(Y),Z,[1 1],'k--','linewidth',2);
            clabel_deterministic(C_temp,ch_temp,'fontsize',10,'LabelSpacing',100);
            drawnow
            clear eff_data;
            for ip = 1:length(ch.EdgePrims)
                eff_data(ip).VertexData = ch.EdgePrims(ip).VertexData;
                eff_data(ip).StripData = ch.EdgePrims(ip).StripData;
            end
            clear eff_ch0
            for ll = 1:length(ch.EdgePrims)
            eff_ch0.EdgePrims(ll).StripData = ch.EdgePrims(ll).StripData;
            eff_ch0.EdgePrims(ll).VertexData = ch.EdgePrims(ll).VertexData;
            end
            
            clabel_deterministic(C,ch,'fontsize',10,'LabelSpacing',200);
            drawnow
            eff_ch = ch;  eff_C = C;
            
            
        else
            clabel_deterministic(C,ch,'fontsize',10,'LabelSpacing',870);
        end
        
    else
        [C,ch] = contour(unique(X),unique(Y),Z,40,'k--','linewidth',1);
    end
    
    
    
    
    if isfield(AR1_AR2_data,metrics{m})
      %  d = plot(AR1_AR2_data.(metrics{m})(:,1),AR1_AR2_data.(metrics{m})(:,2), '.','markersize',4,'markerfacecolor','none','markeredgecolor',repmat(0.4,1,3));
    end
    
    fun = @(x) - Fs{m}(x(1),x(2));
    guesses = [1 0; 4 0.9; 10 0.9;  10 0;  5.5 0.5;  1.5 0];
    [standardized_guesses] = standardize(guesses, limits);
    clear body metric_val
    for g = 1:size(guesses,1)
        %         [body(g,:), metric_val(g)] = fminsearch(fun, standardized_guesses(g,:)',options);
        [body(g,:), metric_val(g)]  = patternsearch(fun, standardized_guesses(g,:)',[],[],[],[],[0 0]',[1 1]',[],options);
    end
    [~,ind] = max(-metric_val);  best_body = unstandardize(body(ind,:),limits);
    
    switch metrics{m}
        case {'eff','construction','tumbling'}
            starh = plot(best_body(1),best_body(2),'kp','markersize',18,'markerfacecolor','k');
        otherwise
            temp = 0.3;
            starh = plot(best_body(1),best_body(2),'p','markersize',18,'markeredgecolor',[temp temp temp],'markerfacecolor',[temp temp temp]);
    end
    
    
    % plot selected optimal shapes
    meshfolder = 'C:\Hull\body meshes\';
    ax = gca;
    
    
    if ismember(m,[1 4 5 6 18])
        switch m
             case 1
                bodies = [10 0];
                positions = [9.75 0.125];  % in data coords
            case 4
                bodies = [10 0.25];
                positions = [9.95 0.29];  % in data coords
            case 5
                bodies = [1.4375 0; 4.5 0.7;];
                positions = [1.6 0.02; 4.6 0.73;];  % in data coords
            case 6
                bodies = [10 0.15];
                positions = [9.8 0.17];
            case 18
                bodies = [1 0;];
                positions = [0.725 0.01;];
        end
        
        for b = 1:size(bodies,1)
            name = ['curved_rod_AR1_',num2str(bodies(b,1)),'_AR2_',num2str(bodies(b,2)),'.dat'];
            [Mesh, Metadata] = load_mesh([meshfolder,name],[],[]);
            if isempty(Mesh)
                continue
            end
            
            [Mesh] = global_inds(Mesh);
            [Mesh] = renumber_Mesh(Mesh);
            
            clear temp_input
            temp_input.performance.nthreads = 8;
            temp_input.accuracy.integration_tol.area.abstol = 1E-6;
            temp_input.accuracy.integration_tol.area.reltol = 1000;
            temp_input.accuracy.integration_tol.area.maxevals = Inf;
            
            temp_input.accuracy.integration_tol.centroid.abstol = 1E-6;
            temp_input.accuracy.integration_tol.centroid.reltol = 1000;
            temp_input.accuracy.integration_tol.centroid.maxevals = Inf;
            
            temp_input.accuracy.integration_tol.volume.abstol = 1E-6;
            temp_input.accuracy.integration_tol.volume.reltol = 1000;
            temp_input.accuracy.integration_tol.volume.maxevals = Inf;
            
            [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
            
            if bodies(b,1) == 1 %sphere degenerate case
                Mesh(1).refpoints = [0 0 0]';   %center
            elseif bodies(b,2) == 0  %not a sphere, but a straight rod degenerate case
                Mesh(1).refpoints = [0 0 0]';   %center
            else  %actually a general curved rod
                Mesh(1).refpoints =  [2*Metadata(1).geom.radius_curv*sin(Metadata(1).geom.nturns / 2 *pi) * cos(Metadata(1).geom.nturns / 2 *pi); ...
                    2*Metadata(1).geom.radius_curv*(sin(Metadata(1).geom.nturns / 2 *pi))^2;
                    0 ];   %center of centerline
            end
            
            
            [Xf, Yf] = ds2nfu(ax,positions(b,1) + 0., positions(b,2) - 0.) ;
            ax_shape = axes('color','none','visible','off','clipping','off','position',[Xf Yf 0.05 0.05]);
            %       set(ax_shape,'visible','on','box','on','xtick','','ytick','');
            axi = axi + 1;    axs(axi) = ax_shape;
            
            
            orig_limits = [1 10; 0 1];
            shape_limits = [-20 20; -20 20];
            orig2shape.x = polyfit(orig_limits(1,:),shape_limits(1,:),1);
            orig2shape.y = polyfit(orig_limits(2,:),shape_limits(2,:),1);
            shape_coord = [polyval(orig2shape.x,positions(b,1)) polyval(orig2shape.y,positions(b,2)) 0]';
            
            [s,e] = plot_mesh(Mesh);
            
            set(e,'edgealpha',0.0);  set(s,'FaceColor',[0.8 0.8 0.8]);
            
            xlim([Mesh(1).Centroid(1) Mesh(1).Centroid(1) + 1]);
            ylim([Mesh(1).Centroid(2) Mesh(1).Centroid(2) + 1]);
            zlim([-1 3.5]);
            
            % hold on
            light
            drawnow
            
        end
        
    end
    
    set(ax,'XTick',[1 2:1:10]);
   
    hold off
    set(ax,'fontsize',14);  fontsizes.labels = 17;
    if ismember(m,[5 4 1])
        xlabel(ax,'SF_1 (elongation)','fontsize',fontsizes.labels);
    end
    if ismember(m,[18 4 1])
        ylabel(ax,'SF_2 (curvature)','fontsize',fontsizes.labels);
    end
    cb.Label.FontSize = fontsizes.labels ;
    %     pause
    %     switch metric
    %         case 'construction'
    %             saveas(handles(m),['C:\Hull\pareto\tasks\',fullnames{m},'.png']);
    %         otherwise
    %             saveas(handles(m),['C:\Hull\pareto\tasks\',displaynames{m},'.png']);
    %     end
    
    if do_save
        switch metrics{m}
            case 'construction'
                saveas(handles(m),['C:\Hull\fitness landscapes\',fullnames{m},'.png']);
                saveas(handles(m),['C:\Hull\fitness landscapes\',fullnames{m},'.fig']);
            otherwise
                saveas(handles(m),['C:\Hull\fitness landscapes\',displaynames{m},'.png']);
                saveas(handles(m),['C:\Hull\fitness landscapes\',displaynames{m},'.fig']);
        end
    end
    
    
    
    
    switch m
        case 5
        
%             for ip = [8 9 11 12]
%                 eff_ch.EdgePrims(ip).VertexData = eff_data(ip).VertexData;
%                 eff_ch.EdgePrims(ip).StripData = eff_data(ip).StripData;
%             end
%             
%             for ip = [13 15 18 19 20]
%                 eff_ch.TextPrims(ip).Visible = 'off';
%             end
drawnow

%             texts = {eff_ch.TextPrims.String};
%             locations = [eff_ch.TextPrims.VertexData];
%            bads =  {'1.0175' , '1.01' ,   '1.015', '1.016',  '1.018', '1.02'};
%            maxnums = [0           1         1         2           1     0    ];
%            for b = 1:length(bads)
%                
%                inds = find(strcmp(bads{b},texts));
%                if length(inds) <= maxnums(b)
%                    continue
%                end
%                sublocations = locations(:,inds);
%                [sorted,inds2] = sort(sublocations(2,:));
%                badlocations = sorted(1:end-maxnums(b));
%                badinds = find(ismember(locations(2,:),badlocations));
%                for bi = 1:length(badinds)
%                    eff_ch.TextPrims(badinds(bi)).Visible = 'off';
%                end
%            end
%            
%            
%             for ll = 1:length(eff_ch.EdgePrims)
%             eff_ch1.EdgePrims(ll).StripData = eff_ch.EdgePrims(ll).StripData;
%             eff_ch1.EdgePrims(ll).VertexData = eff_ch.EdgePrims(ll).VertexData;
%             end
%             
%             linds = [8 9 10 ];
%             for li = 1:length(linds)
%                 eff_ch.EdgePrims(linds(li)).ColorBinding = 'interpolated';
%                 bads = find( eff_ch.EdgePrims(linds(li)).VertexData(1,:) <= 1.6 & eff_ch.EdgePrims(linds(li)).VertexData(2,:) <= 0.436);
%                 eff_ch.EdgePrims(linds(li)).ColorData = uint8(repmat([0 0 0 255]',1,size(eff_ch.EdgePrims(linds(li)).VertexData,2)));
%                 eff_ch.EdgePrims(linds(li)).ColorData(4,bads) = 0;
%             end
               
%                saveas(gcf,'C:\allsync\all papers\curved rod ms\figures\fitness landscapes\efficiency.png');
        case 6
%              saveas(gcf,'C:\allsync\all papers\curved rod ms\figures\fitness landscapes\temporal.png');
        case 18
%              saveas(gcf,'C:\allsync\all papers\curved rod ms\figures\fitness landscapes\construction weighted curv2.png');
    end
    
    set(ax,'position',[  0.132863531353356    0.174766424535667  0.652684567623691   0.741049096194507]);
end


Pareto_data.fullnames = fullnames;
Pareto_data.Fs = Fs;
Pareto_data.metrics = metrics;

return

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
            avg_AR1_AR2 = mean( vertcat(species_data.AR) );  SF2 = avg_AR1_AR2(2);  %SF2 = 0.4
        case 'straight'
            SF2 = 0;
    end
    
    
    temp = standardize( [SF1'  repmat(SF2,length(SF1),1)] , limits);
    
    for i = [1:6 18]
        SF2_transect.(names{n}).(metrics{i}) = Fs{i}(temp);
    end
    
end

%hardcode known optimal SF2 for certain tasks
for i = [1 2 18] %uptake fore-aft construction
    SF2_best.(metrics{i}) = [zeros(size(SF1));  SF2_transect.straight.(metrics{i})'];
end
%% make line graph of curved vs straight performance
figure(300);  clf

names = {'Dm','eff','temporal','uptake','fore_aft','construction'};
% legnames = {'Dispersal','Swimming Efficiency','Chemotaxis (temporal)','Nutrient Uptake','Chemotaxis (spatial)','Construction Ease'};
legnames = {'Dispersal','Swimming Efficiency','Finding Nutrients (temporal strategy)','Nutrient Uptake','Finding Nutrients (spatial strategy)','Construction Cost'};
colors = distinguishable_colors(length(names));
styles = {'-','-','-','--','--','--'};
clear transect_pl
for n = 1:length(names)
    transect_pl(n) = plot(SF2_transect.SF1 , SF2_transect.curved.(names{n}) ./ SF2_transect.straight.(names{n}),styles{n});
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
xlabel('SF_1 (elongation)');  ylabel({'relative performance of','\bf{curved} \rm(SF_2 = 0.15) \bf{vs straight}'});  
set(gca,'fontsize',22);

