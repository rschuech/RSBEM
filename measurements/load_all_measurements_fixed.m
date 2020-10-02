% load dump file, then revert data to data_original, then run later cells -
% no need to run first cell or fix_MicrobeJ_data again


load('C:\allsync\all papers\curved rods\organized\data_dump7.mat')
data = data_original;


% also, at some point H. pylori was manually removed based on crappy fit to
% curved rod shape.  if it still seems to be in data_filtered.good (an outlier above Pareto region around L = 4.5, K = 0.4), do this
% data_filtered.good(120) = [];
%%


from_scratch = false;

% be sure to run this again after all fixing to update size values based on
% fixed length, width

% redo_files = {'Azospirillum irakense sp nov a nitrogen-fixing bacterium associated with rice roots and rhizosphere soil_Page_07_Image_0002_Page_1_Image_0001'};  %the data for these files will be completely overwritten and all fixed data deleted - use if deleting any cells



redo_files = {'Desulfobacter vibrioformis - Lien and Beeder 1997_Page_2_Image_0002_edited',...
    'Desulfonatronum zhilinae sp nov a novel haloalkaliphilic sulfate-reducing bacterium from soda Lake Alginskoe Trans-Baikal Region Russia_Page_4_Image_0001 panel B',...
    'Evaluation of Leptospirillum ferrooxidans for leaching_Page_4_Image_0002 panel A',...
    'Interfacial activity and leaching patterns of Leptospirillum ferrooxidans on pyrite_Page_03_Image_0001',...
    'Leptospirillum ferriphilum - Liu et al 2007_Page_4_Image_0002.txt',...
    'New types of acetate-oxidizing sulfate-reducing Desulfobacter species D hydrogenophilus sp nov D latus sp nov and D curvatus sp nov_Page_3_Image_0002_Page_1_Image_0001 panel D',...
    'vol2C-nnn-910-1114'};

redo_files = {'Leptospirillum ferrooxidans_Giavenoetal2007'};

redo_files = {'Leptospirillum ferrooxidans_Giavenoetal2007' , 'Helicobacter cholecystus - Franklin et al 1996_Page_3_Image_0001_edit_b'};

redo_files = {'Desulfovibrio aespoeensis - Motamedi and Pedersen 1998_Page_2_Image_0002_edited'};

redo_files = {'Caulobacter Segnis - Patel 2014_Page_30_Image_0002'};
%  redo_files = {'Helicobacter pylori - Josenhans et al 1995_Page_07_Image_0002_edited'};
redo_files = {};
% to completely remove all data / mention of an image, need to manually
% find and delete from data struct

folder = 'C:\allsync\all papers\curved rods\organized\all results-14-10\';
master_file = 'C:\allsync\all papers\curved rods\organized\master list checked.xlsx';
SI_file = 'C:\allsync\all papers\curved rod ms\Table S1.xlsx';

[num,txt,raw] = xlsread(master_file,'data');

headers = txt(1,:);
[~,names_col] = ismember('genus species',headers);
[~,files_col] = ismember('exact name of .txt results file',headers);
[~,SF1_col] = ismember('calculated SF1',headers);  [~,SF2_col] = ismember('calculated SF2',headers);
[~,uses_scale_col] = ismember('.txt uses scale?',headers);
[~,multiplication_factor_col] = ismember('.txt multiplication factor for microns',headers);
[~,pixels_per_micron_col] = ismember('pixels per micron conversion factor',headers);
[~,size_info_col] = ismember('size info from text (microns)',headers);


start_col = 'C';  start_row = 3;  %upper left corner of SF1 SF2 rows, cols


files = txt(start_row:end,files_col);  names = txt(start_row:end,names_col);
uses_scale = txt(start_row:end,uses_scale_col);
multiplication_factor = num(1:end,multiplication_factor_col - 2);
pixels_per_micron = num(1:end, pixels_per_micron_col - 2);
length_width_reported = num(1:end, (size_info_col - 2) : (size_info_col - 2 + 5) );  % [min_length  mean_length  max_length  min_width  mean_width  max_width] reported in main text

vars = {'SHAPE.curvature','SHAPE.length','SHAPE.width.mean','FIXED'};  % important variables to save, make sure FIXED is last since it often won't be in the file at all
vars2 = {'curvature','length','width','feret_fixed'};

if from_scratch
    %     clear data
    %     data.file = [];  data.curvature = []; data.length = [];  data.width = [];
    data = [];
    existing_files = {};
else
    existing_files = {data.file};
end


clear cols
SF2write = NaN(length(files),2);  size2write = NaN(length(files),3);  N_cells2write = NaN(length(files),1);

ff = 0;
for f = 1:length(files)
    file = files{f};
    
    fid = fopen([folder,file,'.txt']);
    if fid == -1
        file
        %                 pause
        continue
        
        
    end
    
    
    [~, ff] = ismember(file, existing_files);
    if ff == 0
        ff = length(data) + 1;
    end
    
  
    if ismember(file,redo_files)
        fields = fieldnames(data(ff));
        for fd = 1:length(fields)
            data(ff).(fields{fd}) = [];
        end
    end
    
    % scale shat used to be here
    
    
    temp = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    temp = temp{1}; %always weirdly all in a one-cell cell array
    headers = temp{1};
    headers = textscan(headers,'%s','Delimiter',',');  headers = headers{1};
    [~,var_inds] = ismember(vars,headers); % which columns of the text file correspond to variables of interest
    var_inds(var_inds == 0) = NaN;
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
        
        for dd = 1:length(data_inds)
            
            
            if length(nums) < data_inds(dd) || isnan(data_inds(dd))
                if strcmp(vars2{dd}, 'feret_fixed')
                    data(ff).(vars2{dd})(li) = false;
                else
                    error('variable not found in results file')
                end
            else
                data(ff).(vars2{dd})(li) = nums(data_inds(dd));
            end
            
        end
        
        
    end
    
    
    data(ff).file = file;
    
    
    if contains(uses_scale{f}, 'Yes','IgnoreCase',true)
        if ~isnan(multiplication_factor(f))
            mult_fact = multiplication_factor(f);
        else
            mult_fact = 1;
        end
        if ~isnan(pixels_per_micron(f))
            error('Pixels per micron value found when txt uses scale')
        end
        data(ff).has_scale = true;
    else
        
        
        mult_fact = 1 / pixels_per_micron(f);
        if isnan(mult_fact)
            data(ff).has_scale = false;
        else
            data(ff).has_scale = true;
        end
        
    end
    
    
    
    if isnan(mult_fact)  % if no scale info for this image, rescale to get lengths around 1 just to help fixing code work
        mult_fact = 1 / mean(data(ff).length);
    end
    
    data(ff).length = data(ff).length * mult_fact;
    data(ff).width = data(ff).width * mult_fact;
    data(ff).curvature = data(ff).curvature / mult_fact;
    
    
    if data(ff).has_scale
        if isfield(data,'length_fixed') && length(data(ff).length_fixed) == length(data(ff).length)
            vol = (data(ff).length_fixed - data(ff).width_fixed).*pi.*(data(ff).width_fixed/2).^2 + 4/3*pi*(data(ff).width_fixed/2).^3;
            data(ff).length_scaled = data(ff).length_fixed;  data(ff).width_scaled = data(ff).width_fixed;  data(ff).curvature_scaled = data(ff).curvature_fixed;
        else
            vol = (data(ff).length - data(ff).width).*pi.*(data(ff).width/2).^2 + 4/3*pi*(data(ff).width/2).^3;
            data(ff).length_scaled = data(ff).length;  data(ff).width_scaled = data(ff).width;  data(ff).curvature_scaled = data(ff).curvature;
        end
        
          
        
    else
        
        computed_mean_length = mean(length_width_reported(f,[1 3])); % mean of min, max reported values - would be same as median of two values
        computed_mean_width = mean(length_width_reported(f,[4 6]));  % mean of min, max reported values - would be same as median of two values
        if isnan(length_width_reported(f,2)) % no mean length reported
            mean_length = computed_mean_length;  % use computed mean length based on min, max length
        else
            mean_length = length_width_reported(f,2);
        end
        if isnan(length_width_reported(f,5)) % no mean length reported
            mean_width = computed_mean_width;  % use computed mean length based on min, max length
        else
            mean_width = length_width_reported(f,5);
        end
        
        if isnan(mean_length)   % if only have length or width, use measured SF1 to infer other variable
            % SF1 = length / width
            if isfield(data,'length_fixed') && length(data(ff).length_fixed) == length(data(ff).length)
                mean_length = median(data(ff).length_fixed  ./  data(ff).width_fixed) * mean_width;
            else
                mean_length = median(data(ff).length  ./  data(ff).width) * mean_width;
            end
        elseif isnan(mean_width)
            if isfield(data,'length_fixed') && length(data(ff).length_fixed) == length(data(ff).length)
                mean_width = mean_length ./ median(data(ff).length_fixed  ./  data(ff).width_fixed);
            else
                mean_width = mean_length ./ median(data(ff).length  ./  data(ff).width);
            end
        end
        
        vol = (mean_length - mean_width).*pi.*(mean_width/2).^2 + 4/3*pi*(mean_width/2).^3;
        data(ff).length_scaled = NaN;  data(ff).width_scaled = NaN;  data(ff).curvature_scaled = NaN;
      
        
    end
  
      
        data(ff).sph_dia = ( (vol*3/4/pi).^(1/3) * 2 );
 
    
    
    
    temp = strsplit(names{f});
    
    data(ff).genus = temp{1};
    ind = strfind(temp{2},'.txt');
    if ~isempty(ind)
        temp{2} = temp{2}(1:ind-1);
    end
    data(ff).species = temp{2};
    
    
    
    SF2write(f,:) = [median(data(ff).SF1) median(data(ff).SF2)];
    size2write(f,:) = [median(data(ff).length_scaled)  median(data(ff).width_scaled)  median(data(ff).curvature_scaled)];
    N_cells2write(f) = length(data(ff).SF1);
end  % end files loop


try
%      xlswrite(master_file,SF2write,'data',[start_col,num2str(start_row)]);
%      xlswrite(SI_file,[SF2write size2write  N_cells2write],'data',['D',num2str(4)]);
end



% after loading data, then run fix_MicrobeJ_data.m, then run rest of this
% file
% note that aggregation cell below overwrites data and you need to go back
% to data_original to reload / refix without starting from scratch since
% aggregated data can't really be modified

data_original = data;  % aggregation below will overwrite data

% when loading a dump, must manually set data = data_original and then run
% below cells to regenerate all necessary variables
%% aggregate cells from images that go together (same strain, media, source, etc)
master_file = 'C:\allsync\all papers\curved rods\organized\master list checked.xlsx';

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

data_aggregated = [];  agg_fields = {'curvature','length','width','curvature_fixed','length_fixed','width_fixed','SF1','SF2','sph_dia'};
data_temp = data;  fields = fieldnames(data_temp);  fields2 = ['file','genus','species',agg_fields];  fields2remove = setdiff(fields,fields2);  data_temp = rmfield(data_temp,fields2remove);


files = {data.file}';
for f = 1:length(aggregate.files)  % each group
    [~,inds] = ismember(aggregate.files{f} , files );  %indices of data that match each image to aggregate in current group
    data_aggregated(f).file = {data(inds).file}';
    for a = 1:length(agg_fields)
        data_aggregated(f).(agg_fields{a}) = [data(inds).(agg_fields{a})];  % this is where the aggregation actually happens - the values from multiple data inds are concatenated
    end
    data_aggregated(f).genus = data(inds(1)).genus;  data_aggregated(f).species = data(inds(1)).species;  % presumably all aggregated images within group are same genus, species....
    
end

[non_agg_files,inds] = setdiff(files, agg_files);  % all files in data not in data_aggregated
data_aggregated = [data_temp(inds) data_aggregated];
genus_species = {};
for a = 1:length(data_aggregated)
    genus_species{a} = [data_aggregated(a).genus, ' ' , data_aggregated(a).species];
end
[~,inds] = sort(genus_species);
data_aggregated = data_aggregated(inds);

data = data_aggregated;







genera = unique({data.genus});  all_genera = {data.genus};
species = unique({data.species});  all_species = {data.species};

species_data = [];  genus_data = [];

for g = 1:length(genera)  % each unique genus
    data_genus = data(ismember(all_genera,genera{g})); % all images of this genus
    %     clear means medians samples
    %     for d = 1:length(data_genus) % each image from this genus
    %         SF = [ data_genus(d).SF1' data_genus(d).SF2' ];  % all SF in image
    %         means(d,:) = mean(SF,1);  medians(d,:) = median(SF,1);
    %         samples(d,1) = length(data_genus(d).SF1);
    %
    %     end
    
    %     genus_data(g).genus = genera{g};
    %     genus_data(g).mean_unweighted = mean(means,1);  genus_data(g).median_unweighted = median(medians,1);
    %     genus_data(g).mean_weighted = sum(means.*repmat(samples,1,2),1)/sum(samples);  genus_data(g).median_weighted = sum(medians.*repmat(samples,1,2),1)/sum(samples);
    %     genus_data(g).samples = length(data_genus);  %how many images for this genus (regardless of same or different species)
    %
    %     genus_data(g).SE = std(means,0,1) ./ sqrt(size(means,1));
    %     genus_data(g).SF = means;
    
    species2 = unique({data_genus.species});  %list of all species
    for s = 1:length(species2)  % loop over species
        data_species = data_genus(ismember({data_genus.species}, species2{s}));  %all image data for this species
        clear means medians samples
        for ss = 1:length(data_species)  %images for this species
            SF = [ data_species(ss).SF1' data_species(ss).SF2' ];
            means.SF(ss,:) = mean(SF,1);  medians.SF(ss,:) = median(SF,1);  % means within each (aggregated) image for this species
            means.sph_dia = mean(data_species(ss).sph_dia);  medians.sph_dia = median(data_species(ss).sph_dia);
            samples(ss,1) = length(data_species(ss).SF1);
            
        end
        
        species_data(end+1).genus = genera{g};   species_data(end).species = species2{s};
        vars = {'SF','sph_dia'};
        for v = 1:length(vars)
            
            species_data(end).mean_unweighted.(vars{v}) = mean(means.(vars{v}),1);  
            species_data(end).median_unweighted.(vars{v}) = median(medians.(vars{v}),1);
            species_data(end).mean_weighted.(vars{v}) = sum(means.(vars{v}).*repmat(samples,1,2),1)/sum(samples);  
            species_data(end).median_weighted.(vars{v}) = sum(medians.(vars{v}).*repmat(samples,1,2),1)/sum(samples);  % actually a weighted average of the medians?
        end
        species_data(end).samples = length(data_species);
        species_data(end).n_individuals = length( [data_species.SF1] );  % how many cells in total (over all aggregated images, averaged strains, etc) we have for this species
        
        species_data(end).SE = std(means.SF,0,1) ./ sqrt(size(means.SF,1));
        species_data(end).SF = means.SF;
        species_data(end).sph_dia = means.sph_dia;  % this is janky, should do means or medians all the way down if genus means/medians are ever used
        
        
    end
    
end

genera = unique({species_data.genus});  all_genera = {species_data.genus};
species = unique({species_data.species});  all_species = {species_data.species};
genus_data0 = genus_data;  genus_data = [];

for g = 1:length(genera)
    data_genus = species_data(ismember(all_genera,genera{g}));
    clear means medians samples
    for d = 1:length(data_genus)
        SF = [ data_genus(d).SF(:,1) data_genus(d).SF(:,2) ];
        means.SF(d,:) = mean(SF,1);  medians.SF(d,:) = median(SF,1);
          means.sph_dia = mean(data_genus(ss).sph_dia);  medians.sph_dia = median(data_genus(ss).sph_dia);
        samples(d,1) = size(data_genus(d).SF,1);
        
    end
    
    genus_data(g).genus = genera{g};
    vars = {'SF','sph_dia'};
    for v = 1:length(vars)
        genus_data(g).mean_unweighted.(vars{v}) = mean(means.(vars{v}),1);  genus_data(g).median_unweighted.(vars{v}) = median(medians.(vars{v}),1);
        genus_data(g).mean_weighted.(vars{v}) = sum(means.(vars{v}).*repmat(samples,1,2),1)/sum(samples);  genus_data(g).median_weighted.(vars{v}) = sum(medians.(vars{v}).*repmat(samples,1,2),1)/sum(samples);
        % weighted by number of species in the genus
    end
    genus_data(g).samples = length(data_genus);  %how many species for this genus (regardless of number of images or bugs in each species)
    genus_data(g).n_individuals = sum( [data_genus.n_individuals] );
    genus_data(g).SE = std(means.SF,0,1) ./ sqrt(size(means.SF,1));
    genus_data(g).SF = means.SF;
    
end

% hulls
cutoff_SF1 = true; % cut off outline at SF1 = 10;  NOTE, original data is always used to compute boundary - this only controls whether we cut that boundary at SF1 = 10 for Pareto GOF computation
add_straights = true;  %manually include all straight rods for outline
shrink_factor = 0.5;  %for boundary calculation, 0 is convhull, 1 is as compact as shatlab is willing to do    was 0.75
shrink_factor = 0.5;  %for boundary calculation, 0 is convhull, 1 is as compact as shatlab is willing to do    was 0.75
min_individuals = 1;  % don't count species / genera with fewer than this number of cells in total
% manual_removes = [142];  % Lepto ferro
manual_removes = [];
remove_sp = true;  % remove any sp. species (i.e. known genus, unknown species)
peel_boundary = false;

clear data_filtered
Pareto_data.observed = [];

cases = {'individuals','species','genera'};
% aggregation_method = 'mean_unweighted'; % for species and genera data
aggregation_method = 'median_unweighted'; % for species and genera data

for c = 1:length(cases)
    switch cases{c}
        case 'individuals'
            meas = [ [data.SF1]' [data.SF2]'  ];  % all inidividual data points
        case 'species'
            temp0 = species_data;  temp0(manual_removes) = [];  data_filtered.manually_removed = species_data(manual_removes);
            
            if remove_sp
                species_temp = {temp0.species};  temp0( ismember(species_temp, 'sp.') ) = [];
                species_temp = {species_data.species};   sp_inds = find( ismember(species_temp, 'sp.') );
                
            end
            
            temp = temp0( [temp0.n_individuals] >= min_individuals );   data_meas0 = temp;
            tmp = vertcat(temp.(aggregation_method));  meas0 = vertcat(tmp.SF);  %size0 = vertcat(tmp.sph_dia);
            
            good_inds = setdiff(1:length(species_data), [manual_removes  ,  sp_inds]);
            
            data_filtered.unknown_species = species_data(sp_inds);
            data_filtered.good = species_data(good_inds);
            
            if peel_boundary
                obs_inds = 1:size(meas0,1);
                if add_straights
                    straights = [linspace(1,max(meas0(:,1)),500)' repmat(0,500,1)];
                    straight_inds = (size(meas0,1)+1):((size(meas0,1)) + size(straights,1));
                    meas = [meas0;  straights  ];  % manually insert expected continuous range of straight rods
                else
                    meas = meas0;
                end
                [k,a] = boundary(meas(:,1),meas(:,2),shrink_factor);
                boundary_pt_inds = unique( k( ismember(k, obs_inds) ) );  % inds of meas0 that are observations on the original boundary
                meas = meas0( setdiff(1:size(meas0,1), boundary_pt_inds),:);  % peeled version of meas
                data_filtered.peeled = data_meas0(boundary_pt_inds);
                data_filtered.good = data_meas0( setdiff(1:size(meas0,1), boundary_pt_inds));
            else
                meas = meas0;
                data_filtered.good = data_meas0;
                data_filtered.peeled = [];
            end
            
            
            
        case 'genera'
            temp = genus_data( [genus_data.n_individuals] >= min_individuals );
            tmp = vertcat(temp.(aggregation_method));     meas = vertcat(tmp.SF);
    end
    Pareto_data.observed.(cases{c}).points = meas;
    
    if add_straights
        meas = [meas;   [linspace(1,max(meas(:,1)),500)' repmat(0,500,1)]       ];  % manually insert expected continuous range of straight rods
    end
    [k,a] = boundary(meas(:,1),meas(:,2),shrink_factor);
%     if cutoff_SF1
        rect = [1 1 10 10 1; 0 1 1 0 0];
        [x,y] = polybool('intersection',rect(1,:),rect(2,:),flipud(meas(k,1)),flipud(meas(k,2)));
        Pareto_data.observed.(cases{c}).boundary.cutoff_SF1 = [x' y'];
%     else
        Pareto_data.observed.(cases{c}).boundary.all_data =  [meas(k,1) meas(k,2)];
        ind = find(Pareto_data.observed.(cases{c}).boundary.all_data(:,2) == 0,1,'last');  % index of last point on SF1 axis
        % insert a NaN to break the ordinarily closed boundary outline
        Pareto_data.observed.(cases{c}).boundary.all_data = [Pareto_data.observed.(cases{c}).boundary.all_data(1:ind,:); Pareto_data.observed.(cases{c}).boundary.all_data(ind,1) NaN; Pareto_data.observed.(cases{c}).boundary.all_data(ind+1:end,:) ];
%     end
    
end


temp = vertcat(data_filtered.good.median_unweighted);
temp = vertcat(temp.sph_dia);
figure(491);  histogram(temp);


species_data_simplified = [];
for s = 1:length(species_data)
    species_data_simplified(s).genus = species_data(s).genus;
    species_data_simplified(s).species = species_data(s).species;
    species_data_simplified(s).SF = species_data(s).median_unweighted.SF;
    species_data_simplified(s).sph_dia = species_data(s).median_unweighted.sph_dia;
end
%
% cutoff_SF1 = false;
figure(37);  clf;
if cutoff_SF1
    individual_size = 3;  species_size = 4;
else
    individual_size = 8;  species_size = 6;
end

individuals = plot([data.SF1],[data.SF2],'.','markerfacecolor',repmat(0.4,1,3),'markeredgecolor',repmat(0.4,1,3),'markersize',individual_size,'visible','off');
hold on
temp = vertcat(species_data.(aggregation_method));  temp = vertcat(temp.SF);
species = plot(temp(:,1),temp(:,2),'o','markerfacecolor','b','markersize',species_size);
hold on
% temp = vertcat(genus_data.mean_unweighted);
% genera = plot(temp(:,1),temp(:,2),'ro','markerfacecolor','r','markersize',7);
if ~cutoff_SF1
    bound_individuals = plot(Pareto_data.observed.individuals.boundary(:,1),Pareto_data.observed.individuals.boundary(:,2),'--','linewidth',1,'color',repmat(0.4,1,3),'visible','off');
    bound_species = plot(Pareto_data.observed.species.boundary(:,1),Pareto_data.observed.species.boundary(:,2),'b-','linewidth',1.5);
    % bound_genera = plot(Pareto_data.observed.genera.boundary(:,1),Pareto_data.observed.genera.boundary(:,2),'r-','linewidth',2);
end
if cutoff_SF1
    xlim([1   , 10 + 2E-2]);
else
    xlim([1 27]);
end

ylim([0 - 3E-3 ,  0.8]);
grid;  box off;
xlabel('SF_1 (elongation)'); ylabel('SF_2 (curvature)');

hold off


return
%% plot selected species dots
if cutoff_SF1
    clear selected
    selected.genus =    {'Azospirillum', 'Leptospirillum','Desulfovibrio','Desulfotomaculum', 'Vibrio'      , 'Thioalkalivibrio'  ,     'Helicobacter', 'Desulfovibrio' ,'Chlorobium'    ,  'Ammonifex' ,'Rhodospirillum'  ,     'Vibrio'   , 'Vibrio',           'Heliobacterium'   };
    selected.species =  { 'halopraeferens',   'ferrooxidans', 'africanus',       'acetoxidans' ,     'ruber'          ,  'jannaschii'      ,  'flexispira'  ,    'gigas'    ,'phaeovibrioides' ,'degensii','rubrum','vulnificus' ,  'cholerae',        'undosum'  };
    manual_positions = [    3.049381000642924   0.457093394856654 ; ... Leptospirillum
        2.263249725438614   0.415480024312492 ; ... Helicobacter cholecystus
        1.494623244419264   0.024440079717490 ; ... Desulfotomaculum
        3.259622012941009   0.174287450157401  ; ... Salinivibrio
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
    stopa
    
    clear genus_species SF
    for g = 1:length(data_filtered.good)
        genus_species{g} = [data_filtered.good(g).genus,' ',data_filtered.good(g).species];
        SF(g,:) = data_filtered.good(g).median_unweighted.SF;
    end
    
    figure(37)
    hold on
    psel = plot(SF(:,1),SF(:,2),'bx','MarkerFaceColor','none','MarkerSize',17);
    txtsel = text(SF(:,1)+0.1,SF(:,2),selected.fullname,'fontweight','bold','backgroundcolor','w','margin',eps,'fontsize',12,'fontangle','italic');
    
    
     txtsel = text(SF(:,1)+0.1,SF(:,2),genus_species,'fontweight','bold','backgroundcolor','w','margin',eps,'fontangle','italic','fontsize',11);

end
%%
if cutoff_SF1
    tic;
    nvars = size(Pareto_data.observed.species.points,1);
    
    min_count = 0.5 * nvars; % region must contain at least this fraction of total # species
    
    %  [c,ceq] = bnonlcon(inds,SF1SF2,min_area)
    % [a , k ] = bobj(inds,SF1SF2)
    % shrink_factor = 0.0; % 0 for convhull, 1 for most compact allowable shape
    
    
    % obj = @(inds) bobj(inds,Pareto_data.observed.species.points,shrink_factor);
    % nonlcon = @(inds) bnonlcon(inds,Pareto_data.observed.species.points,shrink_factor, min_count);
    
    % [obj  = ellobj(ell)
    % [c,ceq] = ell_con(ell, SF1SF2,min_count)
    
    nonlcon = @(ell)ell_con(ell, Pareto_data.observed.species.points, min_count);
    options = optimoptions('patternsearch','display','off');
    % guess = [0.2 0.1 5.5 0.05 45*pi/180]';
    lb = [0 0 1 0 0]';  ub = [10 0.5 10 1 2*pi]';
    
    n_guesses = 10000;
    clear ells fvals flags
    % nn = 0;
    ppm = ParforProgressStarter2('ellipse calc', n_guesses);
    parfor n = 1:n_guesses
        guess = lb + (ub - lb).*rand(5,1);
        
        [ells(n,:),fvals(n),flags(n)] = patternsearch(@ellobj,guess,[],[],[],[],lb,ub,nonlcon,options);
        
        %     if flag > 0
        %
        %         nn = nn + 1;
        %         ells(nn,:) = ell_temp;
        %         fvals(nn) = fval_temp;
        %     end
          ppm.increment(n);
    end
    delete(ppm);
    
    ells = ells(flags > 0 , :);  fvals = fvals(flags > 0);
    
    
    [~,inds] = sort(fvals);
    ell = ells(inds(1),:);
    [  ells(inds,:)  fvals(inds)'  ];
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
    
%     median_pt = median( vertcat( species_data.mean_unweighted ) );
  temp = vertcat(data_filtered.good.median_unweighted);   SF = vertcat(temp.SF);  median_pt = median(SF);
    if cutoff_SF1
        mp = plot(median_pt(1),median_pt(2),'gp','markerfacecolor','g','markersize',20);  mp.Color = 'r';
    end
    
    hold off
    toc/60
end
if ~cutoff_SF1
%     legend([individuals species bound_individuals bound_species],{'individuals','species means' ,'individuals boundary','species boundary' },'location','best');
else
    %     legend([individuals species  pl  mp],{'individuals','species means'  , '50% of all species' , 'median species'},'location','best');
    
%     legend([individuals species  pl  mp],{'individual cells','species'  , '50% of all species' , 'average species'},'location','best');
end
% set(gca,'FontSize',18);

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