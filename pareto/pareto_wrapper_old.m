% load E:\Hull\curv_on_1400_grid.mat
% Construction_Ease = 1./mean_abs_curv;
% F_curv = scatteredInterpolant(standardized_X(:), standardized_Y(:), Construction_Ease(:),'natural','none');

goal_1s = {'Construction_Ease','Construction_Ease','Construction_Ease', 'Construction_Ease','Construction_Ease','Construction_Ease','Power_eff','Power_eff',    'Power_eff',         'Power_eff',    'Power_eff'};
goal_2s = {'Power_eff',        'Dm',               'temporal_SN',       'temporal_ability', 'fore_aft_SN',      'fore_aft_ability', 'Dm',       'temporal_SN',  'temporal_ability',  'fore_aft_SN',  'fore_aft_ability'};
%%

styles = {'ks',                    'ks',            'bo',                'bs',                 'ro',               'rs',             'ko',        'ks',           'k^',                'ko',             'ks'};

n_pts = 200;  %how many points along initial guess Pareto line to start searching from


parfor gg = 1:length(goal_1s)
    
    goal_1 = goal_1s{gg};
    goal_2 = goal_2s{gg};
    
    dist_tol = 1E-3;
    max_guesses = 500;
    front = [];
    
    limits = [min(X(:)) max(X(:)); min(Y(:)) max(Y(:))];
    
    
    % c_vec = linspace(0.6, 1, 600);  %for construction ease
    
    
    
    % goal_1 = 'Construction_Ease';  %we step through c_vec contour values for this goal
    % goal_2 = 'fore_aft_ability';   %we try to maximize gaol_2 along each contour of goal_1
    
    
    
    [row,col] = find(Pareto_data.(goal_1).Z == max(Pareto_data.(goal_1).Z(:)));
    AR_1 = [ X(row,col)  Y(row,col) ];
    
    
    [row,col] = find(Pareto_data.(goal_2).Z == max(Pareto_data.(goal_2).Z(:)));
    AR_2 = [ X(row,col)  Y(row,col) ];
    
    
    m =  (AR_2 - AR_1)  ./ ( 1 - 0);  % parameterized slope from t = 0 to 1 from archetype (or max value) for goal 1 to goal 2
    b = AR_1 - m.*0;  % t = 0 at AR_1
    
    
    
    t = linspace(0,1,n_pts)';
    AR = repmat(m,n_pts,1).*repmat(t,1,2) + repmat(b,n_pts,1);
    
    AR_standardized = standardize(AR,limits);
    
    %         nth = 1;
    %     c_vec = Pareto_data.(goal_1).F(standardized_X_subset(1:nth:length(standardized_X_subset)), standardized_Y_subset(1:nth:length(standardized_Y_subset))); %draw a contour through every data point
    
    
    c_vec = Pareto_data.(goal_1).F(AR_standardized(:,1),AR_standardized(:,2));
    
    
    i = 0;
    cc = 0;
    for c = c_vec'
        i = i+1;
        i/length(c_vec)
        
        if isnan(c)
            continue
        end
        
        [ contours] = get_contour(c, unique(X),unique(Y),Pareto_data.(goal_1).Z);  %speed contour
        %     if length(contours.x_c) <= 1
        %         continue
        %     end
        
        if isempty(contours)
            continue
        end
        
        [pp_0] = splining(contours);
        
        if ~isstruct(pp_0)
            continue
        end
        
        x = []; fval = [];
        
        for p = 1:length(pp_0) % in case there are multiple disconnected contours at this contour level
            
            objfun = @(t_0) obj_contours(t_0,  pp_0(p), Pareto_data.(goal_2).F, limits);
            
            wc = 0;
            while true
                wc = wc + 1;
                t_0_guess = rand;
                if ~isnan(objfun(t_0_guess))
                    break
                end
                if wc > max_guesses
                    break
                end
                
            end
            
            if wc > max_guesses
                x(p) = NaN;  fval(p) = NaN;
                continue
            end
            
            [x(p),fval(p)] = patternsearch(objfun,[t_0_guess]',[],[],[],[],[0 ]',[1]' ,optimoptions('patternsearch','display','none','meshtolerance',5E-8,'steptolerance',1E-10,'maxiterations',1200));
            
            if isnan(x)
                stopa
            end
            
        end
        
        %        if wc > max_guesses
        %             continue
        %         end
        
        [fval,ind] = min(fval);
        
        if isnan(fval)
            continue
        end
        
        x = x(ind);    %choose contour that contains best possible F value
        coord = (fnval(pp_0(ind),x))';  %[AR1 AR2] of best point
        
        
        % now check that if we move along other F contour containing
        % current point, we can't find a better value of this F (if we can,
        % this point is not on the Pareto front)
        standardized_coord = standardize(coord,limits);
        
        fval2 = Pareto_data.(goal_2).F(standardized_coord(1),standardized_coord(2));
        [ contours2] = get_contour(fval2, unique(X),unique(Y),Pareto_data.(goal_2).Z);  %speed contour
        %         if length(contours2.x_c) <= 1
        %         continue
        %     end
        if isempty(contours2)
            continue
        end
        
        [pp_0_2 ] = splining(contours2);
        
        x2 = []; fval2 = [];
        for p = 1:length(pp_0_2) % in case there are multiple disconnected contours at this contour level
            objfun2 = @(t_0) obj_contours(t_0,  pp_0_2(p), Pareto_data.(goal_1).F, limits);
            
            
            wc = 0;
            while true
                wc = wc + 1;
                t_0_guess2 = rand;
                if ~isnan(objfun2(t_0_guess2))
                    break
                end
                if wc > max_guesses
                    break
                end
            end
            
            if wc > max_guesses
                x2(p) = NaN;  fval2(p) = NaN;
                continue
            end
            
            [x2(p),fval2(p)] = patternsearch(objfun2,[t_0_guess2]',[],[],[],[],[0 ]',[1 ]' ,optimoptions('patternsearch','display','none','meshtolerance',5E-8,'steptolerance',1E-10,'maxiterations',1200));
        end
        
        [fval2,ind2] = min(fval2);
        
        if isnan(fval2)
            continue
        end
        
        x2 = x2(ind2);    %choose contour that contains best possible F value
        coord2 = (fnval(pp_0_2(ind2),x2))';  %[AR1 AR2] of best point
        
        
        dist = sqrt(sum((coord - coord2).^2));
        
        if dist < dist_tol
            
            cc = cc + 1;
            front(cc,:) = coord;
            
        end
        
    end
    
    fronts(gg).goal_1 = goal_1;  fronts(gg).goal_2 = goal_2;  fronts(gg).points = front;
    
    
end


stopafra
%%
for gg = 1:length(fronts)
    inds = roundn(fronts(gg).points(:,1),-2) == 10   |  fronts(gg).points(:,2) > 0.65;
    fronts(gg).points(inds,:) = [];
end

%%

delete(hans);  clear hans

for gg = finds %%1:length(fronts)
    %     delete(han)
    %     for p = 1:size(fronts(gg).points,1)
    %     han(p) = plot(fronts(gg).points(p,1),fronts(gg).points(p,2),styles{gg},'markerfacecolor',styles{gg}(1),'markersize',8);
    %     drawnow
    %     pause
    %     end
    
    hans(gg) = plot(fronts(gg).points(:,1),fronts(gg).points(:,2),[styles{gg},'-'],'markerfacecolor',styles{gg}(1),'markersize',8);
    drawnow
    fronts(gg).points(1,:)
    pause
end

%%
figure(355)
clear splines legends Hans
finds = [1 2 4 6 7 9 11];
starts = [1 0;  1 0;  1 0;   1 0;   6 0.65;  6 0.65;  6 0.65;];
lasts = [6 0.65;  10 0.3;   10 0.3;  10 0.1;  10 0.56;  10 0.56;  10 0.3];
increasing = {'on','on','on','on','off','off','off'};
decreasing = {'off','off','off','off','on','on','off'};
concaveup = {'off','off','off','off','on','on','off'};
numknots = [9 8 8 6 6 2 2 2];

styles = {'g-','r--','r:','r-','b--','b:','b-'};
try, delete(Hans), end;  try, delete(leg), end;
for f = 1:length(finds)
    points = [starts(f,:); fronts(finds(f)).points; lasts(f,:) ];
    dists = cumsum( [0 ; sqrt( sum( ( points(2:end,:) - points(1:end-1,:) ).^2 , 2) ) ] );
    t = dists / dists(end);
    
    options = slmset('knots',numknots(f),'leftvalue',starts(f,1),'rightvalue',lasts(f,1));
    slm_x = slmengine(t,points(:,1),options);
    options = slmset('knots',numknots(f),'leftvalue',starts(f,2),'rightvalue',lasts(f,2),'minvalue',0,'increasing',increasing{f},'decreasing',decreasing{f},'concaveup',concaveup{f});
    slm_y = slmengine(t,points(:,2),options);
    t_splined = linspace(0,1,200);
    clear splined
    splined(:,1) = slmeval(t_splined,slm_x);   splined(:,2) = slmeval(t_splined,slm_y);
    
    splines{f} = splined;
    %     legends{f}{1} = fronts(finds(f)).goal_1;    legends{f}{2} = fronts(finds(f)).goal_2;
    legends{f} = [  fronts(finds(f)).goal_1,'    ',    fronts(finds(f)).goal_2  ];
    %     delete(han1);  delete(han2)
    %      han1 = plot(fronts(finds(f)).points(:,1),fronts(finds(f)).points(:,2),[styles{finds(f)},'-'],'markerfacecolor',styles{finds(f)}(1),'markersize',8);
    %      han2 = plot(splined(:,1),splined(:,2),'k-','linewidth',2);
    %      drawnow
    %      pause
    
    Hans(f) = plot(splined(:,1),splined(:,2),styles{f},'linewidth',2);  hold on;
    
    %
end
xlim([1 10 + 1.5E-2]);  ylim([0 1]);  grid on;  %just enough xlim to show data point at AR1 = 10

try, delete(ch), end;  [C,ch] = contour(unique(X),unique(Y),Pareto_data.Construction_Ease.Z,[0.8 0.8],'-.','color',repmat(0.7,1,3));
set(ch,'linewidth',2);

set(Hans([1 4 7]),'linewidth',5);  %Construction, Efficiency, Fore-Aft appear to be the overall Pareto region

try, delete(leg); end
% leg = legend([Hans, ch],{'Construction Ease - Swimming Efficiency','Construction Ease - Dispersal','Construction Ease - Temporal Chemotaxis', ...
%     'Construction Ease - Fore-Aft Chemotaxis','Swimming Efficiency - Dispersal','Swimming Efficiency - Temporal Chemotaxis',...
%     'Swimming Efficiency - Fore-Aft Chemotaxis', '80% of best Construction Ease'},'interpreter','none');
% leg = legend([Hans],{'Construction Ease - Swimming Efficiency','Construction Ease - Dispersal','Construction Ease - Temporal Chemotaxis', ...
%     'Construction Ease - Fore-Aft Chemotaxis','Swimming Efficiency - Dispersal','Swimming Efficiency - Temporal Chemotaxis',...
%     'Swimming Efficiency - Fore-Aft Chemotaxis'},'interpreter','none');
% set(leg,'fontsize',fontsize-6);


leg = legendflex([Hans],{'Construction Ease - Swimming Efficiency','Construction Ease - Dispersal','Construction Ease - Temporal Chemotaxis', ...
    'Construction Ease - Fore-Aft Chemotaxis','Swimming Efficiency - Dispersal','Swimming Efficiency - Temporal Chemotaxis',...
    'Swimming Efficiency - Fore-Aft Chemotaxis'},'ncol',2,'fontsize',fontsize - 6);

set(leg,'Position',[ 528       691.23          895        115.3]);
set(gca,'fontsize',18);
fontsize = 22;
xlabel('AR_1 (elongation)','fontsize',fontsize);   ylabel('AR_2 (curvature)','fontsize',fontsize);

try, delete(starh), end
starh = plot([1 6],[0 0.65],'kp','markerfacecolor','k','markersize',24);

set(gca,'xtick',1:10);
patpts = [ splines{1}; splines{7}; flipud(splines{4}) ];
try, delete(pat), end;
pat = patch(patpts(:,1),patpts(:,2),'k','facealpha',0.15,'edgecolor','none');

try, delete(hullh), end;
xtemp = X(~isnan(Pareto_data.Power_eff.Z));  ytemp = Y(~isnan(Pareto_data.Power_eff.Z));
k = convhull(xtemp,ytemp,'simplify',true);
hullpts = [xtemp(k) ytemp(k)];
hullpts(hullpts(:,1) == 10 & hullpts(:,2) <= 0.95 - 1E-2 , :) = []; % remove right vertical boundary
hullpts(hullpts(:,2) == 0 & hullpts(:,1) > 1 , :) = [];  %remove bottom horizontal boundary
hullpts(1,:) = [];  %get rid of annoying connecting line
hullh = plot(hullpts(:,1),hullpts(:,2),':','linewidth',2,'color',[repmat(0.7,1,3)]);

%%
print('-dpng','-r600',['E:\Hull\select_dumps\','Pareto_headliner','.png']);


%%

load('E:\Hull\CFD03\dino_code\pareto\swimming_efficiency.mat')
%%
c = 0.001;



[ contours] = get_contour(c, unique(X),unique(Y),Z);

[pp_0, ppder_0 ] = splining(contours);
%%

load('E:\Hull\CFD03\dino_code\pareto\Dm.mat')

%%
% c = 0.99
% t_0 = 0.6475;
% c_1 = 0.3202;
% t_1 = 0.42;

t_0 = 0.;
c_1 = 0.5;
t_1 = 1;

% curvature
V = 1;
[mean_curv, max_curv, mean_abs_curv, mean_max_curv] = curved_rod_curvature(X, Y, V);

% [standardized_X, standardized_Y] = standardize(X,Y,[min(X(:)) max(X(:)); min(Y(:)) max(Y(:))]);  %already have this shat from swimming_efficiency.m

F_curv = scatteredInterpolant(standardized_X(:), standardized_Y(:), 1./(mean_abs_curv(:)),'natural','none');


X2 = X;  Y2 = Y;  %use same X and Y as for swimming eff
Z2 = 1./mean_abs_curv;  Z2 = Z2 / max(Z2(:)); %high curvature is bad, low curvature is good

x_1 = unique(X2);  y_1 = unique(Y2);  Z_1 = Z2;

% [obj] = obj_contours_wrapper(t_0, c_1, t_1,  pp_0, ppder_0, x_1, y_1, Z_1)
%%

objfun = @(t_0) obj_contours(t_0,  pp_0, F_curv);

t_0_guess = 0.5;

x = patternsearch(objfun,[t_0_guess]',[],[],[],[],[0 ]',[1 ]' ,optimoptions('patternsearch','display','iter','meshtolerance',5E-8,'steptolerance',1E-10,'maxiterations',1200,'initialmeshsize',4))


objfun = @(x) obj_contours_wrapper(x(1), x(2), x(3),  pp_0, ppder_0, x_1, y_1, Z_1);

x = patternsearch(objfun,[t_0 c_1 t_1]',[],[],[],[],[0 0 0]',[1 1 1]' ,optimoptions('patternsearch','display','iter','meshtolerance',5E-8,'steptolerance',1E-10,'maxiterations',1200,'initialmeshsize',4))

%%

opts = optimoptions(@fmincon,'Algorithm','interior-point','display','iter');
problem = createOptimProblem('fmincon','objective',...
    objfun,'x0',[t_0 c_1 t_1]','lb',[0 0 0]','ub',[1 1 1]','options',opts);
gs = GlobalSearch;
[x,f] = run(gs,problem)


%%
t_sampled = linspace(0,1,10000);
contour_0 = fnval(pp_0,t_sampled);
start_pt = fnval(pp_0,0); last_pt = fnval(pp_0,1);

pt_0 = fnval(pp_0,x(1));  der_0 = fnval(ppder_0, x(1));  %der_0 = der_0(2) / der_0(1);

[ contours_1] = get_contour(x(2), x_1,y_1,Z_1);
[pp_1, ppder_1 ] = splining(contours_1);
contour_1 = fnval(pp_1,t_sampled);
pt_1 = fnval(pp_1,x(3));  der_1 = fnval(ppder_1, x(3));  %der_1 = der_1(2) / der_1(1);

[pt_0 pt_1]
[der_0/sqrt(sum(der_0.^2)) der_1/sqrt(sum(der_1.^2))]

figure(69)
try, delete(poop), end
poop = plot(contour_0(1,:),contour_0(2,:),'r-',pt_0(1),pt_0(2),'ro',   contour_1(1,:),contour_1(2,:),'b-',pt_1(1),pt_1(2),'bo');
set(poop,'linewidth',2);
xlim([1 12]);  ylim([0 1]);


%%
nAR1 = 150;  nAR2 = 150;
%  nAR1 = 60;  nAR2 = 60;
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
          fontsizes.labels = 18;
        fontsizes.axes = 18;
        fontsizes.title = 15;
        fontsizes.caxis = 26;
    fontsizes.legend = 22;
    
    
goals = {'Construction_Ease','Power_eff','Dm','temporal_SN','fore_aft_SN', 'tumbling'  ,  'uptake'};
goals_titles = {'Construction Ease','Swimming Efficiency','Dispersal','Temporal S/N','Fore-Aft S/N','Tumbling Ease'  ,  'Uptake'};
names_titles = {'Construction','Swimming','Dispersal','Temporal','Fore-Aft','Tumbling' , 'Uptake'};

% goals = {'translation_ease','rotation_resistance'};
% goals_titles = {'translation_ease','rotation_resistance'};
% names_tiles = {'translation_ease','rotation_resistance'};

goals = {'Construction_Ease','Power_eff','fore_aft_SN'};
goals_titles = {'Construction Ease','Swimming Efficiency','Fore-Aft S/N'};
names_titles = {'Construction','Swimming','Fore-Aft'};


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
    goals = goal_combos{go};
    
%     ranks = []; perfs = [];
%     for g  = 1:length(goals)
%         ranks = cat(3,ranks,Pareto_data.(goals{g}).rank);
%         perfs(:,g) = interped.(goals{g})(:);
%         
%     end
    

%     siz = size(Pareto_data.Power_eff.rank);
% %     optimal = ones(size(Pareto_data.Power_eff.rank));  % 1 for optimal, 0 for suboptimal, NaN for outside parameter space
%       optimal = zeros(size(Pareto_data.Power_eff.rank));  
%     optimal( any(isnan(ranks),3) ) = NaN;
%     
%     nl = size(ranks,1)*size(ranks,2);
%     tic;
%     for i = 1:nl
% %         i / nl
%         if isnan(optimal(i))
%             continue
%         end
%         [ri,ci] = ind2sub(siz,i);
%         
% %         if ~( ri == 6 && ci == 31)
% %             continue
% %         else
% %             pausefuck
% %         end
%         
%         for ii = 1:nl
%             
% %             if optimal(ii) == 0 || isnan(optimal(ii))
% %                 continue
% %             end
%               if isnan(optimal(ii))
%                 continue
%             end
%             
%             
%             [rii,cii] = ind2sub(siz,ii);
%             
%             tol = 0.005;
% %        is_beaten =    any( (perfs(ii,:) - perfs(i,:)) ./ ( (perfs(ii,:) + perfs(i,:))/2 ) >= tol );
% %        beats = any( (perfs(i,:) - perfs(ii,:)) ./ ( (perfs(ii,:) + perfs(i,:))/2 ) >= tol );
%        
%            is_beaten =    any( (perfs(ii,:) - perfs(i,:)) ./ ( (perfs(ii,:) + perfs(i,:))/2 ) >= tol );
%        beats = any( (perfs(i,:) - perfs(ii,:)) ./ ( (perfs(ii,:) + perfs(i,:))/2 ) >= tol );
%        
% %        if is_beaten && ~beats
% %             optimal(i) = false;
% %                 break;
% %        end
%        
%          if ~is_beaten && beats
%             optimal(i) = true;
%                 break;
%        end
%             
%             
%             % for g = 1:length(goals)
% %             if all( (ranks(ri,ci,:)) < (ranks(rii,cii,:)) )  %another point exists that's better in all goals than this point
% % %                 if ri == 8 && ci == 30
% % %                     fuckdo
% % %                 end
% %                 optimal(i) = false;
% %                 break;
% %             end
%         end
%     end
%     toc
    
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


   figure(825)
    clf
       
   %  [Contour_data,buckling_handle] = contour(unique(X),unique(Y),Pareto_data.buckling_force_mean.Z,repmat(1.1,1,2),'b-','linewidth',1.5);
%  [Contour_data2,buckling_handle2] = contour(unique(X),unique(Y),Pareto_data.buckling_force_max.Z,[1.25 1.3 1.35 1.4 ],'r-','linewidth',1.5);
 %[Contour_data2,buckling_handle2] = contour(unique(X),unique(Y),Pareto_data.buckling_force_max.Z,repmat(1.3,1,2),'b-','linewidth',1.5);
% max of 1.3 is nearly identical to mean of 1.1, but 1.3 sounds better...
buckling_poly = [Contour_data2(:,2:end) [10; 1]  [Contour_data2(1,2); 1] [Contour_data2(:,2)] ];
% [Contour_data,ease_handle] = contour(unique(X),unique(Y),Pareto_data.Construction_Ease.Z,[0.8 0.75],'k--','linewidth',2);
buckling_too_high  = find(inpolygon(XPc,YPc,buckling_poly(1,:),buckling_poly(2,:)));

% [Contour_data,contour_handle] = contour(unique(X),unique(Y),Pareto_data.fore_aft_SN.Z,[1 1],'k--','linewidth',1);
Contour_data = contourcs(unique(X),unique(Y),Pareto_data.fore_aft_SN.Z,[1 1]);
% contour_handle = plot(Contour_data(2).X,Contour_data(2).Y,'k--','linewidth',1);


hold on
   

    [ p, optimal_inds] = paretoFront(rounded );
    
    % remove points with buckling over the assumed hard limit from the optimal
% set

%  optimal_inds = setdiff(optimal_inds,buckling_too_high);
    
    
    Optimal = NaN(size(XP));
    for oi = 1:length(optimal_inds)
        Optimal( XP == XPc(optimal_inds(oi)) & YP == YPc(optimal_inds(oi)) ) = 1;
    end
   
    A_pareto = dA * length(optimal_inds);  % approx area of Pareto region
%     ninds = find(any(isnan(perfs),2));
%     inds = setdiff(inds,ninds);
%     po = plot(XPc(inds),YPc(inds),'bo','markerfacecolor','b','markersize',8); xlim([1 10]); ylim([0 1]);
%    in_data_hull = inpolygon(XPc,YPc,observed.individuals(:,1),observed.individuals(:,2)); %returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
     in_data_hull = inpolygon(XPc,YPc,observed.species(:,1),observed.species(:,2)); %returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    
%    po = plot(XPc(optimal == 1),YPc(optimal == 1),'bo','markerfacecolor','b','markersize',8); xlim([1 10]); ylim([0 1]);
    
 
    pc = pcolor(XP,YP,(Optimal));  
 
    set(pc,'FaceAlpha',0.15);
   colormap(gray(20));
    shading flat
   hold on

%    temp  = surf2patch(pc1);
%    pc = patch(temp);


    % plot optimal pts
%      try, delete(po);, end;  po = plot(X_Y_unq(optimal==1,1),X_Y_unq(optimal==1,2),'bo','markerfacecolor','b','markersize',8);
        in_optimal = intersect( find(in_data_hull) , optimal_inds );
        out_optimal = intersect( find(~in_data_hull) , optimal_inds);
        in_suboptimal = setdiff( find(in_data_hull) , optimal_inds );

%       poo = plot(XPc(in_suboptimal),YPc(in_suboptimal),'ko','markersize',10);
      hold on
%        po = plot(XPc(out_optimal),YPc(out_optimal),'bo','markerfacecolor','b','markersize',8);
   temp = vertcat(plot_data.mean_unweighted);
all_meas = plot(meas(:,1),meas(:,2),'ko','markerfacecolor','k','markersize',1);  % individual dots
% bound_individuals = plot(observed.individuals(:,1),observed.individuals(:,2),'k--','linewidth',1.5); %individuals hull
  species_meas = plot(temp(:,1),temp(:,2),'ko','markerfacecolor','k','markersize',10); %species circles
 
  %bound_species = plot(observed.species(:,1),observed.species(:,2),'k-','linewidth',1.5); %species hull
  bound_species = plot(observed.species(1:2,1),observed.species(1:2,2),'k-','linewidth',1.5); %species hull
  plot(observed.species(3:end,1),observed.species(3:end,2),'k-','linewidth',1.5)

%     delete(pat)
%     shape = alphaShape(standardize(X_Y_unq(optimal==1,:), limits),0.255);
%     shape.Points = unstandardize(shape.Points,limits);
%     [bf,P] = boundaryFacets(shape);
%      pat = patch('faces',bf,'vertices',P);

C = dA*length(in_optimal);
B = dA*length(out_optimal);
A = dA*length(in_suboptimal);

GOF = C/(B+C) * C/(A+C);

GOFs(go) = GOF;

  meas = [ [data.AR1]' [data.AR2]'  ];  % all inidividual data points
edges = {  1:0.5:26  ,  0:0.05:1   };
[N,C] = hist3(meas,'edges',edges);
N(N == 0) = NaN;

N2 = padarray(N,[1 1],0,'post');
C2{1} = (C{1}(1)-0.5/2):0.5:(C{1}(end)+0.5/2);
C2{2} = (C{2}(1)-0.05/2):0.05:(C{2}(end)+0.05/2);

%  pcolor(C2{1},C2{2},N2'); xlim([1 10]); ylim([0 1]);  shading flat;  colorbar;  hold on
%  cblabel('# individual cells');



   
% temp = X_Y_unq(optimal==1,:);
%     temp2 = standardize(temp, limits);
%     k = boundary(temp2(:,1),temp2(:,2),0.75);
%     
%     bound = plot(temp(k,1),temp(k,2),'k-','linewidth',2);
   
   
    axis normal

    % plot all body data pts
%     d = plot(X_Y_unq(subset,1),X_Y_unq(subset,2), '.','markersize',3,'markerfacecolor','none','markeredgecolor','k');
    hold on
    





    % plot edges of parameter space
    try, delete(hullh), end;
    xtemp = X(~isnan(Pareto_data.Power_eff.Z));  ytemp = Y(~isnan(Pareto_data.Power_eff.Z));
    k = convhull(xtemp,ytemp,'simplify',true);
    hullpts = [xtemp(k) ytemp(k)];
    hullpts(hullpts(:,1) == 10 & hullpts(:,2) <= 0.95 - 1E-2 , :) = []; % remove right vertical boundary
    hullpts(hullpts(:,2) == 0 & hullpts(:,1) > 1 , :) = [];  %remove bottom horizontal boundary
    hullpts(1,:) = [];  %get rid of annoying connecting line
    hullh = plot(hullpts(:,1),hullpts(:,2),':','linewidth',2,'color',[repmat(0.7,1,3)]);
    
    
    
    [~,ind] = near(Contour_data(2).X(1),hullpts(:,1));
    verts = [flipud(hullpts(1:ind,:)); [Contour_data(3).X' Contour_data(3).Y']; flipud([Contour_data(2).X' Contour_data(2).Y'])];
    faces = 1:length(verts);
    infeasible = patch('Faces',faces,'Vertices',verts,'visible','off');
%     infeasible_h = hatchfill2(infeasible,'HatchSpacing',10,'HatchLineWidth',0.5);
    infeasible_h = hatchfill2(infeasible,'cross');
    set(infeasible_h,'visible','on','color',repmat(0.8,1,3));
    
    xlim([1 10+1.5E-2]);
    ylim([0 1]);
    grid on
        set(gca,'fontsize',fontsizes.axes);
        
%     title({'Goals Considered:',title_combos{go}},'fontsize',fontsizes.title,'interpreter','none');
    title({title_combos{go} , ['Goodness of Fit = ',num2str(GOF)]},'fontsize',13);

    
    
    xlabel('AR_1 (elongation)','fontsize',fontsizes.labels)
    ylabel('AR_2 (curvature)','fontsize',fontsizes.labels)
    set(gcf,'Position',[ 680          82        1051         896]);
    
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
set(legend_h,'Position',[153.16 677.3 266 132.5])
%     hold off
%    pause
    drawnow
%     if GOF > 0.6
% print('-dpng','-r300',['E:\Hull\pareto\current\','GOF ',num2str(GOF*100,'%4.0f'),'%   ',name_combos{go},'.png']);  %  curved tail
%     end
      %pause
     % shateroonie
%  print('-dpng','-r300',['E:\Hull\graphs2\',name_combos{go},' straight opt tail','.png']);  %  curved tail
 
% print('-dpng','-r300',['E:\Hull\pareto\','GOF ',num2str(GOF*100,'%4.0f'),'%   ',name_combos{go},'.png']);  %  curved tail
   
end