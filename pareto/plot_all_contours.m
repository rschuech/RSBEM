clear options;  options.Spacing = [X(1,2)-X(1,1), Y(2,1)-Y(1,1)];
smoothed = smoothn(Z,30, options);
%%
for cl = [0.1:0.015:0.95 0.96:0.005:0.99 0.992:0.001:0.999 ]
    [temp] = contourcs(unique(X),unique(Y),smoothed,[cl cl]);
    lengths = [temp.Length];
    [~,ind] = max(lengths);
    if ~isempty(ind)
    plot(temp(ind).X,temp(ind).Y,'k-');
    end
end

%%
for cl = [0.1:0.01:1 ]
    [temp] = contourcs(unique(X2),unique(Y2),Z2,[cl cl]);
    lengths = [temp.Length];
    [~,ind] = max(lengths);
    if ~isempty(ind)
    plot(temp(ind).X,temp(ind).Y,'k--');
    end
end



%%

[mean_curv, max_curv, mean_abs_curv, mean_max_curv] = curved_rod_curvature(X, Y, 1);
%%
figure(45);  
	p = panel();
	p.pack(2, 2);
    
    p(1, 1).select();
temp = 1./mean_curv;  pcolor(X,Y,temp/max(temp(:)));  shading interp; xlim([1 12]); ylim([0 1]); caxis([0 1]);  hold on;  [C,h] = contour(unique(X),unique(Y),temp,75,'k--');
title('1/(surface averaged mean curvature)')  
p(1, 2).select();
temp = 1./mean_abs_curv;  pcolor(X,Y,temp/max(temp(:)));  shading interp;  xlim([1 12]); ylim([0 1]);  caxis([0 1]);  hold on;  [C,h] = contour(unique(X),unique(Y),temp,75,'k--');
title('1/(surface averaged mean absolute curvature)')  
p(2, 1).select();
temp = 1./mean_max_curv;  pcolor(X,Y,temp/max(temp(:)));  shading interp;  xlim([1 12]); ylim([0 1]);  caxis([0 1]);  hold on;  [C,h] = contour(unique(X),unique(Y),temp,40,'k--');
title('1/(surface averaged max curvature)')  
p(2, 2).select();
temp = 1./max_curv;  pcolor(X,Y,temp/max(temp(:)));  shading interp;  xlim([1 12]); ylim([0 1]);  caxis([0 1]);  hold on;  [C,h] = contour(unique(X),unique(Y),temp,40,'k--');
title('1/(global max curvature)')  

p.fontsize = 12;


% p.de.margin = 1;

%%

figure(83)

temp = 1./mean_abs_curv;  pcolor(X,Y,temp/max(temp(:)));  shading interp;  xlim([1 12]); ylim([0 1]);  caxis([0 1]);  hold on;  [C,h] = contour(unique(X),unique(Y),temp,75,'k--');
title('1/(surface averaged mean absolute curvature) and smoothed power efficiency') 

clear options;  options.Spacing = [X(1,2)-X(1,1), Y(2,1)-Y(1,1)];
smoothed = smoothn(Z,100, options);
[C,h] = contour(unique(X),unique(Y),smoothed,300,'k-');