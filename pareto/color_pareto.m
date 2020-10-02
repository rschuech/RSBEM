


Optimal2 = false(size(Optimal));

Optimal2(Optimal == 1) = true;

Optimal_color = NaN([size(Optimal) 3]);

clear opt_metric_normalized


for task = 1:3
    
    metric = squeeze(Pareto_data.interped(:,:,goals(task)));
    
    
    opt_metric = metric(Optimal2);
    
    metric_range = [min(opt_metric(:)) max(opt_metric(:))];
    
    opt_metric_normalized = ( opt_metric - metric_range(1) ) / (metric_range(2) - metric_range(1));
    full_normalized = NaN(size(Optimal));
    full_normalized(Optimal2) = opt_metric_normalized;
    
    
    Optimal_color(:,:,task) = full_normalized;
    
    
    
end


%%
figure(384);  clf;
xy=[0 0; 1 0; 0.5 sqrt(3)/2];
col=[0 0 1; 0 1 0; 1 0 0];
kp = patch('Vertices',xy, 'Faces',[1:size(xy,1)], 'EdgeColor','none','FaceVertexCData', col,'FaceColor','interp');
axis off
text(xy(1,1), xy(1,2)-0.05,'Construction','HorizontalAlignment', 'center')
text(xy(2,1), xy(2,2)-0.05,'Chemotaxis','HorizontalAlignment', 'center')
text(xy(3,1), xy(3,2)+0.05,'Efficiency','HorizontalAlignment', 'center')
kp.FaceAlpha = 0.75;

