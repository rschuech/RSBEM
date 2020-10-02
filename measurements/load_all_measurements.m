

folder = 'C:\allsync\all papers\curved rods\organized\all results-14-10\';
master_file = 'C:\allsync\all papers\curved rods\organized\master list checked.xlsx';
[num,txt,raw] = xlsread(master_file,'data');

headers = txt(1,:);
[~,names_col] = ismember('genus species',headers);
[~,files_col] = ismember('exact name of .txt results file',headers);
[~,SF1_col] = ismember('calculated SF1',headers);  [~,SF2_col] = ismember('calculated SF2',headers);
[~,uses_scale_col] = ismember('.txt uses scale?',headers);
[~,multiplication_factor_col] = ismember('.txt multiplication factor for microns',headers);
[~,pixels_per_micron_col] = ismember('pixels per micron conversion factor',headers);


start_col = 'C';  start_row = 3;  %upper left corner of SF1 SF2 rows, cols


files = txt(start_row:end,files_col);  names = txt(start_row:end,names_col);

vars = {'SHAPE.curvature','SHAPE.length','SHAPE.width.mean','FIXED'};  % important variables to save, make sure FIXED is last since it often won't be in the file at all


clear data cols
SF2write = NaN(length(files),2);

ff = 0;
for f = 1:length(files)
    file = files{f}; 
    
    fid = fopen([folder,file,'.txt']);
    if fid == -1
                file
%                 pause
        continue
   
      
    end
    temp = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    temp = temp{1}; %always weirdly all in a one-cell cell array
    headers = temp{1};
    headers = textscan(headers,'%s','Delimiter',',');  headers = headers{1};
    [~,var_inds] = ismember(vars,headers); % which columns of the text file correspond to variables of interest
    [~,position_ind] = ismember('POSITION',headers);
    data_inds = var_inds - position_ind + 1;  % convert indices of headers to indices of the extracted subset of values from each data line (which doesn't include all variables in the headers)
    
    temp = temp(2:end); %remove header row
    data_temp = [];
    for li = 1:length(temp)
        line = temp{li};
        ind = strfind(line,')');
        subline = line(ind(end)+2:end);
        nums = textscan(subline,'%f','Delimiter',',');
        nums = nums{1};
        
        data_temp(:,li) = nums(var_inds);
    end
    
    ff = ff + 1;
    data(ff).file = file;
    data(ff).curvature = data_temp(4,:);
    data(ff).length = data_temp(5,:);
    data(ff).width = data_temp(6,:);
    
    
    data(ff).AR1 = data(ff).length ./ data(ff).width;
    radius_curv = 1./ data(ff).curvature;
    data(ff).AR2 = data(ff).length ./ (2*pi*radius_curv);
    
    temp = strsplit(names{f});
    
    data(ff).genus = temp{1};
    ind = strfind(temp{2},'.txt');
    if ~isempty(ind)
        temp{2} = temp{2}(1:ind-1);
    end
    data(ff).species = temp{2};
    
    
    
    SF2write(f,:) = [mean(data(ff).AR1) mean(data(ff).AR2)];
    
    
end
try
xlswrite(master_file,SF2write,'data',[start_col,num2str(start_row)]);
end
%% aggregate cells from images that go together (same strain, media, source, etc)
[num,txt,raw] = xlsread(master_file,'aggregate');
headers = txt(1,:);
[~,names_col] = ismember('genus species',headers);
[~,files_col] = ismember('name of .txt results files to be combined',headers);
[~,strains_col] = ismember('strain (if multiple aggregated strains)',headers);

start_row = 2;
agg_files = txt(start_row:end,files_col);  names = txt(start_row:end,names_col);
strains = txt(start_row:end,strains_col);

clear aggregate
aggregate.indices = {1};  aggregate.names = {names{1}};  aggregate.files{1} = {agg_files{1}};
last_name = names{1};  last_strain = strains{1}; %assume first row is from first group
for n = 2:(length(names))
    if strcmpi( names{n} , last_name ) && strcmpi( strains{n} , last_strain ) % current name, strain matches prev row
        aggregate.indices{end}(end+1) = n;  aggregate.files{end}{end+1,1} = agg_files{n};
    else
        aggregate.indices{end+1} = n;  aggregate.files{end+1}{1} = agg_files{n}; aggregate.names{end+1} = names{n};
        last_name = names{n};  last_strain = strains{n};
    end
end

data_aggregated = [];  agg_fields = {'curvature','length','width','AR1','AR2'};

files = {data.file}';
for f = 1:length(aggregate.files)  % each group
    [~,inds] = ismember(aggregate.files{f} , files );  %indices of data that match each image to aggregate in current group
    data_aggregated(f).file = {data(inds).file}';
    for a = 1:length(agg_fields)
        data_aggregated(f).(agg_fields{a}) = [data(inds).(agg_fields{a})];
    end
    data_aggregated(f).genus = data(inds(1)).genus;  data_aggregated(f).species = data(inds(1)).species;  % presumably all aggregated images within group are same genus, species....
    
end

[non_agg_files,inds] = setdiff(files, agg_files);  % all files in data not in data_aggregated
data_aggregated = [data(inds) data_aggregated];
genus_species = {};
for a = 1:length(data_aggregated)
    genus_species{a} = [data_aggregated(a).genus, ' ' , data_aggregated(a).species];
end
[~,inds] = sort(genus_species);
data_aggregated = data_aggregated(inds);

% data_original = data;  data = data_aggregated;
%%

genera = unique({data.genus});  all_genera = {data.genus};
species = unique({data.species});  all_species = {data.species};

species_data = [];  genus_data = [];

for g = 1:length(genera)  % each unique genus
    data_genus = data(ismember(all_genera,genera{g})); % all images of this genus
%     clear means medians samples
%     for d = 1:length(data_genus) % each image from this genus
%         AR = [ data_genus(d).AR1' data_genus(d).AR2' ];  % all AR in image
%         means(d,:) = mean(AR,1);  medians(d,:) = median(AR,1);
%         samples(d,1) = length(data_genus(d).AR1);
%         
%     end
    
%     genus_data(g).genus = genera{g};
%     genus_data(g).mean_unweighted = mean(means,1);  genus_data(g).median_unweighted = median(medians,1);
%     genus_data(g).mean_weighted = sum(means.*repmat(samples,1,2),1)/sum(samples);  genus_data(g).median_weighted = sum(medians.*repmat(samples,1,2),1)/sum(samples);
%     genus_data(g).samples = length(data_genus);  %how many images for this genus (regardless of same or different species)
%     
%     genus_data(g).SE = std(means,0,1) ./ sqrt(size(means,1));
%     genus_data(g).AR = means;
    
    species2 = unique({data_genus.species});  %list of all species
    for s = 1:length(species2)  % loop over species
        data_species = data_genus(ismember({data_genus.species}, species2{s}));  %all image data for this species
        clear means medians samples
        for ss = 1:length(data_species)  %images for this species
            AR = [ data_species(ss).AR1' data_species(ss).AR2' ];
            means(ss,:) = mean(AR,1);  medians(ss,:) = median(AR,1);  % means for each image for this species
            samples(ss,1) = length(data_species(ss).AR1);
            
        end
        
        species_data(end+1).genus = genera{g};   species_data(end).species = species2{s};
        species_data(end).mean_unweighted = mean(means,1);  species_data(end).median_unweighted = median(medians,1);
        species_data(end).mean_weighted = sum(means.*repmat(samples,1,2),1)/sum(samples);  species_data(end).median_weighted = sum(medians.*repmat(samples,1,2),1)/sum(samples);
        species_data(end).samples = length(data_species);
        
        species_data(end).SE = std(means,0,1) ./ sqrt(size(means,1));
        species_data(end).AR = means;
        
        
    end
    
end

genera = unique({species_data.genus});  all_genera = {species_data.genus};
species = unique({species_data.species});  all_species = {species_data.species};
genus_data0 = genus_data;  genus_data = [];

for g = 1:length(genera)
    data_genus = species_data(ismember(all_genera,genera{g}));
    clear means medians samples
    for d = 1:length(data_genus)
        AR = [ data_genus(d).AR(:,1) data_genus(d).AR(:,2) ];
        means(d,:) = mean(AR,1);  medians(d,:) = median(AR,1);
        samples(d,1) = size(data_genus(d).AR,1);
        
    end
    
    genus_data(g).genus = genera{g};
    genus_data(g).mean_unweighted = mean(means,1);  genus_data(g).median_unweighted = median(medians,1);
    genus_data(g).mean_weighted = sum(means.*repmat(samples,1,2),1)/sum(samples);  genus_data(g).median_weighted = sum(medians.*repmat(samples,1,2),1)/sum(samples);
    % weighted by number of species in the genus
    genus_data(g).samples = length(data_genus);  %how many species for this genus (regardless of number of images or bugs in each species)
    
    genus_data(g).SE = std(means,0,1) ./ sqrt(size(means,1));
    genus_data(g).AR = means;
    
end

%% hulls
cutoff_AR1 = true; % cut off data and outline at AR1 = 10
add_straights = true;  %manually include all straight rods for outline
shrink_factor = 0.5;  %for boundary calculation, 0 is convhull, 1 is as compact as shatlab is willing to do    was 0.75

cases = {'individuals','species','genera'};
aggregation_method = 'mean_unweighted'; % for species and genera data

for c = 1:length(cases)
    switch cases{c}
        case 'individuals'
            meas = [ [data.AR1]' [data.AR2]'  ];  % all inidividual data points
        case 'species'
            meas = vertcat(species_data.(aggregation_method));
        case 'genera'
            meas = vertcat(genus_data.(aggregation_method));
    end
    Pareto_data.observed.(cases{c}).points = meas;
    
    if add_straights
        meas = [meas;   [linspace(1,max(meas(:,1)),500)' repmat(0,500,1)]       ];  % manually insert expected continuous range of straight rods
    end
    [k,a] = boundary(meas(:,1),meas(:,2),shrink_factor);
    if cutoff_AR1
        rect = [1 1 10 10 1; 0 1 1 0 0];
        [x,y] = polybool('intersection',rect(1,:),rect(2,:),flipud(meas(k,1)),flipud(meas(k,2)));
        Pareto_data.observed.(cases{c}).boundary = [x' y'];
    else
        Pareto_data.observed.(cases{c}).boundary =  [meas(k,1) meas(k,2)];
        ind = find(Pareto_data.observed.(cases{c}).boundary(:,2) == 0,1,'last');  % index of last point on SF1 axis
        % insert a NaN to break the ordinarily closed boundary outline
        Pareto_data.observed.(cases{c}).boundary = [Pareto_data.observed.(cases{c}).boundary(1:ind,:); Pareto_data.observed.(cases{c}).boundary(ind,1) NaN; Pareto_data.observed.(cases{c}).boundary(ind+1:end,:) ];
    end
    
end





%%
figure(37);  clf;
if cutoff_AR1
    individual_size = 3;  species_size = 4;
else
    individual_size = 8;  species_size = 6;
end
    
individuals = plot([data.AR1],[data.AR2],'.','markerfacecolor',repmat(0.4,1,3),'markeredgecolor',repmat(0.4,1,3),'markersize',individual_size);
hold on
temp = vertcat(species_data.mean_unweighted);
species = plot(temp(:,1),temp(:,2),'o','markerfacecolor','b','markersize',species_size);
hold on
% temp = vertcat(genus_data.mean_unweighted);
% genera = plot(temp(:,1),temp(:,2),'ro','markerfacecolor','r','markersize',7);
if ~cutoff_AR1
bound_individuals = plot(Pareto_data.observed.individuals.boundary(:,1),Pareto_data.observed.individuals.boundary(:,2),'--','linewidth',1,'color',repmat(0.4,1,3));
bound_species = plot(Pareto_data.observed.species.boundary(:,1),Pareto_data.observed.species.boundary(:,2),'b-','linewidth',1.5);
% bound_genera = plot(Pareto_data.observed.genera.boundary(:,1),Pareto_data.observed.genera.boundary(:,2),'r-','linewidth',2);
end
if cutoff_AR1
    xlim([1   , 10 + 2E-2]);
else
    xlim([1 27]);
end

ylim([0 - 3E-3 ,  0.55]);
grid;  box off;
xlabel('SF_1 (elongation)'); ylabel('SF_2 (curvature)');

hold off

%%
do_histogram = false;
if do_histogram
hold on
try delete(hi); delete(p); end;  [hi,hc] = honeycomb(species.XData,species.YData,round([40 7]*.65));  set(hi,'FaceAlpha',0.3,'edgecolor','none'); %1.15 sort of worked (non unique solution)
f = ones(sum(hc.isIncluded),1);
A = - hc.counts(hc.isIncluded)' ;
b = - 0.5 * sum(hc.counts);  
lb = zeros(size(f));  ub = ones(size(lb));
intcon = 1:length(f);

% region = logical( intlinprog(f, intcon, A, b, [],[],lb,ub) );

min_count = 0.5 * sum(hc.counts);  


subx = hc.XVertices(:,hc.isIncluded);  suby = hc.YVertices(:,hc.isIncluded);
subcounts = hc.counts(hc.isIncluded);

[sorted,inds] = sort(subcounts,'descend');

sorted_cumsum = cumsum(sorted);
last_ind = find(sorted_cumsum >= min_count,1,'first');
    if sorted(last_ind+1) == sorted(last_ind)
        disp('Non-unique solution');
    end
    
    region = false(size(subcounts));
    region( inds(1:last_ind) ) = true;


p = patch(subx(:, (region)),suby(:, (region)),'k', 'edgecolor','k','facecolor','none','linewidth',2,'linestyle','dashed');

sum(subcounts(region)) / sum(subcounts) * 100
figure(36); 

% sort( subcounts(region),'descend')'
% sort(subcounts,'descend')'
end
% return
%% plot selected species dots
if cutoff_AR1
clear selected
selected.genus =    {'Bdellovibrio' ,'Leptospirillum','Helicobacter','Desulfotomaculum', 'Salinivibrio'    ,    'Nitrosovibrio'  , 'Thioalkalivibrio'  ,     'Helicobacter', 'Desulfovibrio' ,'Chlorobium'    ,  'Ammonifex' ,'Pseudomonas'  ,     'Vibrio'   , 'Vibrio',           'Heliobacterium'   };
selected.species =  {'sp.'  ,   'ferrooxidans',  'cholecystus',       'acetoxidans' ,     'costicola' ,              'sp.'            ,  'jannaschii'      ,  'flexispira'  ,    'gigas'    ,'phaeovibrioides' ,'degensii','pseudoalcaligenes','vulnificus' ,  'cholerae',        'undosum'  };
manual_positions = [    2.244645489442075   0.252669782263099   ; ... Bdellovibrio
                       3.049381000642924   0.457093394856654 ; ... Leptospirillum
                        2.263249725438614   0.415480024312492 ; ... Helicobacter cholecystus
                        1.494623244419264   0.024440079717490 ; ... Desulfotomaculum
                        3.259622012941009   0.174287450157401  ; ... Salinivibrio
                        5.641004392974755   0.358961650818647  ; ... Nitrosovibrio
                         9.257225551029260   0.149963371487689 ; ... Thioalkalivibrio
                         9.453578286985008   0.050826233245763 ; ... Helicobacter flexispira
                         7.669561549148290   0.229375681435339 ; ... Desulfovibrio
                         1.037650589973532   0.184188792485720  ; ... Chlorobium
                         5.293886743687068   0.025769131047570  ; ... Ammonifex
                         3.178110917929055   0.042256139192282 ; ... Pseudomonas
                         3.939088677855212   0.277338075216555 ; ... Vibrio vulnificus
                         5.576917306510497   0.222852456210746  ; ... Vibrio cholerae
                         7.473431414971873   0.087851682603609 ; ... Heliobacterium
                         ];
                         

for ss = 1:length(selected.genus)
    for s = 1:length(species_data)
        if strcmp( [species_data(s).genus,' ',species_data(s).species] , [selected.genus{ss},' ',selected.species{ss}])
            selected.species_ind(ss) = s;
            break
        end
    end
    selected.data_inds{ss} = [];
    for s = 1:length(data)
        if strcmp( [data(s).genus,' ',data(s).species] , [selected.genus{ss},' ',selected.species{ss}]);
            selected.data_inds{ss}(end+1) = s;
        end
    end
    selected.fullname{ss} = [selected.genus{ss} , ' ',selected.species{ss}];
end
temp = vertcat(species_data(selected.species_ind).mean_unweighted);
figure(37)
hold on
psel = plot(temp(:,1),temp(:,2),'bo','MarkerFaceColor','b','MarkerSize',7);
% txtsel = text(temp(:,1)+0.1,temp(:,2),selected.fullname,'fontweight','bold','backgroundcolor','w','margin',4);
txtsel = text(manual_positions(:,1),manual_positions(:,2),selected.fullname,'fontweight','bold','backgroundcolor','w','margin',eps,'fontangle','italic','fontsize',11);
end
%%
if cutoff_AR1
tic;
nvars = size(Pareto_data.observed.species.points,1);

min_count = 0.5 * nvars; % region must contain at least this fraction of total # species

%  [c,ceq] = bnonlcon(inds,AR1AR2,min_area)
% [a , k ] = bobj(inds,AR1AR2)
% shrink_factor = 0.0; % 0 for convhull, 1 for most compact allowable shape


% obj = @(inds) bobj(inds,Pareto_data.observed.species.points,shrink_factor);
% nonlcon = @(inds) bnonlcon(inds,Pareto_data.observed.species.points,shrink_factor, min_count);

% [obj  = ellobj(ell)
% [c,ceq] = ell_con(ell, AR1AR2,min_count)

nonlcon = @(ell)ell_con(ell, Pareto_data.observed.species.points, min_count);
options = optimoptions('patternsearch','display','off');
% guess = [0.2 0.1 5.5 0.05 45*pi/180]';
lb = [0 0 1 0 0]';  ub = [10 0.5 10 1 2*pi]'; 

n_guesses = 3000;  
clear ells fvals flags
% nn = 0;
parfor n = 1:n_guesses
    guess = lb + (ub - lb).*rand(5,1);
    
    [ells(n,:),fvals(n),flags(n)] = patternsearch(@ellobj,guess,[],[],[],[],lb,ub,nonlcon,options);
    
%     if flag > 0
%         
%         nn = nn + 1;
%         ells(nn,:) = ell_temp;
%         fvals(nn) = fval_temp;
%     end
    
end

ells = ells(flags > 0 , :);  fvals = fvals(flags > 0); 


[~,inds] = sort(fvals);
ell = ells(inds(1),:);
[  ells(inds,:)  fvals(inds)'  ]
% options = optimoptions('ga','display','iter','PopulationSize',400,'UseParallel',true);
% binds = ga(obj,nvars,[],[],[],[],zeros(1,nvars),ones(1,nvars),nonlcon,1:nvars,options);

% [obj_val, k] = obj(binds);
% min_count - nonlcon(ell) 
[~,~,is_inside] = nonlcon(ell);
[min_count sum(is_inside)]
% subpts = Pareto_data.observed.species.points(logical(binds),:);
% pl = plot(subpts(k,1),subpts(k,2),'b-','linewidth',2);
% https://math.stackexchange.com/questions/941490/whats-the-parametric-equation-for-the-general-form-of-an-ellipse-rotated-by-any
t = linspace(0,1,200);
x = ell(3) +    ell(1)* cos(2*pi*t) *cos(ell(5)) - ell(2)* sin(2*pi*t) *sin(ell(5));
y = ell(4) +    ell(1)* cos(2*pi*t) *sin(ell(5)) + ell(2) *sin(2*pi*t) *cos(ell(5));
% figure(543)
hold on
pl = plot(x,y,'r-','linewidth',1);

median_pt = median( vertcat( species_data.mean_unweighted ) );
if cutoff_AR1
    mp = plot(median_pt(1),median_pt(2),'gp','markerfacecolor','g','markersize',20);  mp.Color = 'r';
end

hold off
toc/60
end
if ~cutoff_AR1
 legend([individuals species bound_individuals bound_species],{'individuals','species means' ,'individuals boundary','species boundary' },'location','best');
else
%     legend([individuals species  pl  mp],{'individuals','species means'  , '50% of all species' , 'median species'},'location','best');

    legend([individuals species  pl  mp],{'individual cells','species'  , '50% of all species' , 'average species'},'location','best');
end
set(gca,'FontSize',18);

return
%%
for i = 1:10
    ell2 = ells(inds(i),:);
t = linspace(0,1,200);
x = ell2(3) +    ell2(1)* cos(2*pi*t) *cos(ell2(5)) - ell2(2)* sin(2*pi*t) *sin(ell2(5));
y = ell2(4) +    ell2(1)* cos(2*pi*t) *sin(ell2(5)) + ell2(2) *sin(2*pi*t) *cos(ell2(5));
% figure(543)
hold on
pl = plot(x,y,'k-','linewidth',2);

end