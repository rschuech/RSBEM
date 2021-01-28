function [colors] = colordata(n,mapname,limits,data)


%n = 20;  % # unique colors; higher --> smoother colormap

% eval(['colors = colormap(',mapname,'(n));']); %will pick from these colors to color each data point

colors = whitejet(n);


colorlo = limits(1);  %data value to be mapped to first color.  data below this saturates
colorhi = limits(2); %data value to be mapped to last color     data above this saturates

%datalo = -20;  %min value of random test data
%datahi = 150;   %max value of random test data

%data = datalo + (datahi-datalo)*rand(20,1);  %random data between datalo and datahi for testing purposes

map = linspace(colorlo,colorhi,n);  %linearly spaced map from data vals to color vals, with colorbar limits colorlo and colorhi

colorinds = interp1(map,1:n,data,'nearest','extrap');  %find nearest matching color index for each data point.  then use this nearest color to color the data.  
%note:  cannot use "linear" since the goal is to output integer indices.  Instead, to get smoother colorbar, increase n

colors = colors(colorinds,:);
%generally, colorinds matrix is matched to data vector -  use colorinds
%to pick the color for each data point
% 
% close
% for i = 1:length(data)
% plot(i,data(i),'o','markerfacecolor',colors(colorinds(i),:),'markeredgecolor',colors(colorinds(i),:))
% hold on
% end
% hold off
% 
% colorbar
% 
% caxis([colorlo colorhi])
% colormap(colors)