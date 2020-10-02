



fraction = 1;  %what fraction of orig data range to extend both sides of plotted interpolant by

subplot(1,3,1)
% buffer = (max(plotdata.amp) - min(plotdata.amp)) * fraction;
% amp = linspace(min(plotdata.amp)-buffer,max(plotdata.amp)+buffer,400);

amp = linspace(min([Results.amp]),max([Results.amp]),400);

% xx = [ repmat(body(1),1,length(amp)); repmat(body(2),1,length(amp)); amp; repmat(lambda_best,1,length(amp)); repmat(nlambda_best,1,length(amp));];
% ef = net(xx) / max(eff.curved_rod.power_effs);
hold on
%pl = plot(amp,ef,'b--');

% vline(net_guess(1),'r--');
% xx = [ repmat(body(1),1,length(amp)); repmat(body(2),1,length(amp)); amp; repmat(net_guess(2),1,length(amp)); repmat(net_guess(3),1,length(amp));];
% ef = net(xx) / max([Results.power_effs]);
% pl = plot(amp,ef,'r--');

vline(pm_guess(1),'r-');
xx = [ repmat(body(1),1,length(amp)); repmat(body(2),1,length(amp)); amp; repmat(pm_guess(2),1,length(amp)); repmat(pm_guess(3),1,length(amp));];
ef = polyvaln(pm, xx') / max([Results.Power_eff]);
pl = plot(amp,ef,'r-');


hold off



subplot(1,3,2)
% buffer = (max(plotdata.lambda) - min(plotdata.lambda)) * fraction;
% lambda = linspace(min(plotdata.lambda)-buffer,max(plotdata.lambda)+buffer,400);

lambda = linspace(min([Results.lambda]),max([Results.lambda]),400);


% xx = [ repmat(body(1),1,length(lambda)); repmat(body(2),1,length(lambda)); repmat(amp_best,1,length(lambda)); lambda; repmat(nlambda_best,1,length(lambda));];
% ef = net(xx) / max(eff.curved_rod.power_effs);
hold on
% pl = plot(lambda,ef,'b--');
% vline(net_guess(2),'r--');
% xx = [ repmat(body(1),1,length(lambda)); repmat(body(2),1,length(lambda)); repmat(net_guess(1),1,length(lambda)); lambda; repmat(net_guess(3),1,length(lambda));];
% ef = net(xx) / max([Results.Power_eff]);
% pl = plot(lambda,ef,'r--');

vline(pm_guess(2),'r-');
xx = [ repmat(body(1),1,length(lambda)); repmat(body(2),1,length(lambda)); repmat(pm_guess(1),1,length(lambda)); lambda; repmat(pm_guess(3),1,length(lambda));];
ef = polyvaln(pm, xx') / max([Results.Power_eff]);
pl = plot(lambda,ef,'r-');

hold off

subplot(1,3,3)
% buffer = (max(plotdata.nlambda) - min(plotdata.nlambda)) * fraction;
% nlambda = linspace(min(plotdata.nlambda)-buffer,max(plotdata.nlambda)+buffer,400);

nlambda = linspace(min([Results.nlambda]),max([Results.nlambda]),400);

% xx = [ repmat(body(1),1,length(nlambda)); repmat(body(2),1,length(nlambda)); repmat(amp_best,1,length(nlambda)); repmat(lambda_best,1,length(nlambda)); nlambda;];
% ef = net(xx) / max(eff.curved_rod.power_effs);
hold on
% pl = plot(nlambda,ef,'b--');
% vline(net_guess(3),'r--');
% xx = [ repmat(body(1),1,length(nlambda)); repmat(body(2),1,length(nlambda)); repmat(net_guess(1),1,length(nlambda)); repmat(net_guess(2),1,length(nlambda)); nlambda;];
% ef = net(xx) / max([Results.Power_eff]);
% pl = plot(nlambda,ef,'r--');

vline(pm_guess(3),'r-');
xx = [ repmat(body(1),1,length(nlambda)); repmat(body(2),1,length(nlambda)); repmat(pm_guess(1),1,length(nlambda)); repmat(pm_guess(2),1,length(nlambda)); nlambda;];
ef = polyvaln(pm, xx') / max([Results.Power_eff]);
pl = plot(nlambda,ef,'r-');

hold off