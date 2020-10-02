
folder = 'E:\allsync\all papers\curved rods\images\microbej\';

files = dir([folder,'*.txt']);
files = {files.name};

clear data cols


for f = 1:length(files)
    file = files{f};
    temp = csvimport([folder,file]);
    
    headers = temp(1,:);
    
    %     cols.AR = find(strcmp('SHAPE.aspectRatio',headers));   %this is for a fitted ellipse, so useless
    cols.curvature = find(strcmp('SHAPE.curvature',headers));
    cols.length = find(strcmp('SHAPE.length',headers));
    cols.width = find(strcmp('SHAPE.width.mean',headers));
    
    
    data(f).file = file;
    fields = fieldnames(cols);
    for fi = 1:length(fields)
        
        data(f).(fields{fi}) = [temp{2:end,cols.(fields{fi})}];
    end
    
    data(f).AR1 = data(f).length ./ data(f).width;
    radius_curv = 1./ data(f).curvature;
    data(f).AR2 = data(f).length ./ (2*pi*radius_curv);
    
    temp = strsplit(file);
    
    data(f).genus = temp{1};
    ind = strfind(temp{2},'.txt');
    if ~isempty(ind)
        temp{2} = temp{2}(1:ind-1);
    end
    data(f).species = temp{2};
end

data(end+1).file = 'Caulobacter crescentus.txt';
data(end).length = 3.25;
data(end).curvature = 1/4.44;
data(end).width = 0.74;
data(end).AR1 = data(end).length ./ data(end).width;
radius_curv = 1./ data(end).curvature;
data(end).AR2 = data(end).length ./ (2*pi*radius_curv);
data(end).genus = 'Caulobacter';
data(end).species = 'crescentus';

%%

genera = unique({data.genus});  all_genera = {data.genus};
species = unique({data.species});  all_species = {data.species};

species_data = [];  genus_data = [];

for g = 1:length(genera)
    data_genus = data(ismember(all_genera,genera{g}));
    clear means medians samples
    for d = 1:length(data_genus)
        AR = [ data_genus(d).AR1' data_genus(d).AR2' ];
        means(d,:) = mean(AR,1);  medians(d,:) = median(AR,1);
        samples(d,1) = length(data_genus(d).AR1);
        
    end
    
    genus_data(g).genus = genera{g};
    genus_data(g).mean_unweighted = mean(means,1);  genus_data(g).median_unweighted = median(medians,1);
    genus_data(g).mean_weighted = sum(means.*repmat(samples,1,2),1)/sum(samples);  genus_data(g).median_weighted = sum(medians.*repmat(samples,1,2),1)/sum(samples);
    genus_data(g).samples = length(data_genus);  %how many images for this genus (regardless of same or different species)
    
      genus_data(g).SE = std(means,0,1) ./ sqrt(size(means,1));
    genus_data(g).AR = means;
    
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
%%
genus_list = {species_data.genus};  species_list = {species_data.species};
[~,temp] = xlsread('E:\allsync\all papers\curved rods\species_flagellation.xls');
for i = 1:size(temp,1)
    ind = find( strcmp(temp{i,1},genus_list) & strcmp(temp{i,2},species_list) );
    species_data(ind).flagellation = temp{i,3};
end


%%
flagellation = {species_data.flagellation}';
filter = strcmp('polar (single)',flagellation);
plot_data = species_data(filter);

% plot_data = species_data;

return
% plot_data = data;
figure(344)
try, delete(hmeas), end;  clear hmeas
clear han
for i = 1:length(plot_data)
     hmeas(i) = plot((plot_data(i).AR1),(plot_data(i).AR2),'ko','markerfacecolor','k','markersize',2);
%       hmeas(i) = plot((plot_data(i).mean_unweighted(1)),(plot_data(i).mean_unweighted(2)),'ko','markerfacecolor','k','markersize',5);

    
    %         if mean(data(i).AR1) > 7
    %             data(i).file
    %             i
    %         end
    
end

%%

 meas = [ [data.AR1]' [data.AR2]'  ]; 
 genus_meas = vertcat(genus_data.mean_unweighted)

figure(394)
subplot(2,1,1)
hist(log10(meas(:,1)),50);  xlim([0 2]); grid on
title('individual cells')
subplot(2,1,2)
hist(log10(genus_meas(:,1)));  xlim([0 2]); grid on
title('genus averages (unweighted)')
xlabel('log(AR_1)')

