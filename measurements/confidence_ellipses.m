[x, y] = meshgrid(-0.5:1:2.5);
% Note I had to pad C to get the 3-by-3 grid you were looking for
C = [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0];
figure
pcolor(x, y, C)


%%

figure;

for i = 1:length(data)
  
        if length(data(i).AR1) == 1
            %plot(data(i).AR1,data(i).AR2,'bo','MarkerFaceColor','none');
        elseif length(data(i).AR1) == 2
           % plot(data(i).AR1,data(i).AR2,'ro','MarkerFaceColor','none');
        elseif length(data(i).AR1) >= 3 &&  length(data(i).AR1) <= 15
             plot(data(i).AR1,data(i).AR2,'ko','MarkerFaceColor','none');
        else
            
            temp = error_ellipse(cov(data(i).AR1,data(i).AR2),[mean(data(i).AR1) mean(data(i).AR2)],'conf',0.95);
            set(temp,'color','k');
    end
    hold on
end

xlim([1 Inf]); ylim([0 1]);  xlabel('AR_1'); ylabel('AR_2');

%%

plot_data = genus_data;
figure(24);

for i = 1:length(plot_data)
  
        if size(plot_data(i).AR,1) == 1
            plot(plot_data(i).AR(1),plot_data(i).AR(2),'bo','MarkerFaceColor','b');
%         elseif length(data(i).AR1) == 2
%            % plot(data(i).AR1,data(i).AR2,'ro','MarkerFaceColor','none');
        elseif size(plot_data(i).AR,1) >= 2 %&&  length(data(i).AR1) <= 15000
%              plot(data(i).AR1,data(i).AR2,'ko','MarkerFaceColor','none');
       % else
       tstar = -tinv((1-0.95)/2,plot_data(i).samples - 1);
           pe =   ploterr(plot_data(i).AR(:,1),plot_data(i).AR(:,2),tstar*plot_data(i).SE(1),tstar*plot_data(i).SE(2)) ;
           set(pe(1),'linestyle','none','marker','o','markerfacecolor','k','color','k');
           set(pe(2),'color','k');  set(pe(3),'color','k');
             %plots X vs. Y with x error bars [X-XE X+XE] and y
%   error bars [Y-YE Y+YE].

%             temp = error_ellipse(cov(plot_data(i).AR(:,1),plot_data(i).AR(:,2)),plot_data(i).mean_unweighted,'conf',0.95);
%             set(temp,'color','k');
    end
    hold on
end

xlim([1 10]); ylim([0 1]);  xlabel('AR_1'); ylabel('AR_2');