%
%     Average over each image and use per-image value as your individual measurement.
%     Plot mean values for individual species, with vertical and horizontal 95% CI estimates as whiskers.
%     If only one species, plot that species, but as a simple point.
%
%     Later we might want to size-code the points based on sample size or another parameter.
%     Ignore strains etc. for now.
%
%     [Repeat on separate figure for Genera rather than species.]
%
%
%
%         also, alpha hull of entire individual bug point cloud


% try, delete(hmeas); end;
%     for i = 1:length(plot_data)
%         %      hmeas(i) = plot((plot_data(i).AR1),(plot_data(i).AR2),'bo','markerfacecolor','b','markersize',1);
%         hmeas(i) = plot((plot_data(i).mean_unweighted(1)),(plot_data(i).mean_unweighted(2)),'ko','markersize',6);
%     end

%     figure(12)



%%


meas = [ [data.AR1]' [data.AR2]'  ];  % all inidividual data points
% hmeas = plot(meas(:,1),meas(:,2),'ko','markerfacecolor','k','markersize',1);
%%
hull = convhull(meas(:,1),meas(:,2));
% hh = plot(meas(hull,1),meas(hull,2),'k:');
%%
% plot_data = species_data;

cutoff_AR1 = true; % cut off data and outline at AR1 = 10
add_straights = true;  %manually include all straight up to AR1 = 10 for outline

n = 390;
figure(n)
clf
edges = {  1:0.5:26  ,  0:0.05:1   };
[N,C] = hist3(meas,edges);
N(N == 0) = NaN;
pcolor(C{1},C{2},N');
if cutoff_AR1
    xlim([1 10]);
else
    xlim([1 Inf]);
end
ylim([0 1]);  shading interp;  colorbar;  hold on

meas = [ [data.AR1]' [data.AR2]'  ];  % all inidividual data points
all_meas = plot(meas(:,1),meas(:,2),'ko','markerfacecolor','k','markersize',1);
if add_straights
    meas = [meas;   [linspace(1,10,20)' repmat(0,20,1)]       ];  % manually insert expected continuous range of straight rods
end
k = boundary(meas(:,1),meas(:,2),0.75);
if cutoff_AR1
    rect = [1 1 10 10 1; 0 1 1 0 0];
    [x,y] = polybool('intersection',rect(1,:),rect(2,:),flipud(meas(k,1)),flipud(meas(k,2)));
    observed.individuals = [x' y'];
else
    observed.individuals = [meas(k,1) meas(k,2)];
end

bound_individuals = plot(observed.individuals(:,1),observed.individuals(:,2),'k--','linewidth',2);

temp = vertcat(genus_data.mean_unweighted);
% species_meas = plot(temp(:,1),temp(:,2),'ko','markersize',6);
% set(species_meas,'MarkerEdgeColor',[repmat(0,1,3)]);
% set(species_meas,'MarkerFaceColor',[repmat(0,1,3)]);

if add_straights
    temp = [temp;   [linspace(1,10,20)' repmat(0,20,1)]       ];  % manually insert expected continuous range of straight rods
end
k = boundary(temp(:,1),temp(:,2),0.75);

if cutoff_AR1
    [x,y] = polybool('intersection',rect(1,:),rect(2,:),flipud(temp(k,1)),flipud(temp(k,2)));
    
    observed.species = [x' y'];
else
    observed.species = [temp(k,1) temp(k,2)];
end
bound_species = plot(observed.species(:,1),observed.species(:,2),'k-','linewidth',2);


figure(n)

