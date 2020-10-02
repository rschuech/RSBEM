
start_angle = 0;
end_angle = pi/4;
range = end_angle - start_angle;


phases = [0              pi/2                  pi        3/4*2*pi              2*pi];
angles = [start_angle    range*(1/8)        range/2     (1-1/8)*range         end_angle];

% options = slmset('interiorknots','free','knots',3,'leftvalue',start_angle,'rightvalue',end_angle,'leftslope',0,'rightslope',0);
options = slmset('leftvalue',start_angle,'rightvalue',end_angle,'leftslope',0,'rightslope',0,'concaveup',[0 pi],'concavedown',[pi 2*pi]);
flipping_spline = slmengine(phases,angles,options);

figure(22)
plot(phases,angles,'o');
hold on
phases2 = linspace(min(phases),max(phases),2000);
plot(phases2,slmeval(phases2,flipping_spline),'-');
hold off
grid on


save('C:\Users\rudi\Desktop\RD\dino turning\tail_flipping_spline.mat','flipping_spline');


%%

angles = fliplr(angles);  

% options = slmset('interiorknots','free','knots',3,'leftvalue',start_angle,'rightvalue',end_angle,'leftslope',0,'rightslope',0);
options = slmset('leftvalue',end_angle,'rightvalue',start_angle,'leftslope',0,'rightslope',0,'concaveup',[pi 2*pi],'concavedown',[0 pi]);
unflipping_spline = slmengine(phases,angles,options);

figure(23)
plot(phases,angles,'o');
hold on
phases2 = linspace(min(phases),max(phases),2000);
plot(phases2,slmeval(phases2,unflipping_spline),'-');
hold off
grid on


save('C:\Users\rudi\Desktop\RD\dino turning\tail_unflipping_spline.mat','unflipping_spline');
