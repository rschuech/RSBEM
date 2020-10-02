
% skips = [1.2 0.375; 1.75 0.075; 1.75 0.125; 1.975 0.625; 2.5 0.425; 2.5 0.475; 2.5 0.525; 2.5 0.575; 2.5 0.625; 2.5 0.675; 2.5 0.725; 2.5 0.775; ...
%     3 0.825; 3.5 0.175; 3.5 0.425; 3.5 0.475; 3.5 0.525; 3.5 0.575; 3.5 0.625; 3.5 0.675; 3.5 0.725; 3.5 0.775; 3.5 0.825; 3.5 0.8775; ...
%     4 0.55; 4 0.575; 4 0.625; 4 0.86; 4 0.96; 4.5 0.275; 4.5 0.425; 4.5 0.475; 4.5 0.525; 4.5 0.725; 4.5 0.775; 4.5 0.825; 4.5 0.875; ...
%     4.5 0.925; 5 0.86; 5 0.92; 5.5 0.275; 5.5 0.425; 5.5 0.475; 5.5 0.525; 5.5 0.725; 5.5 0.775; 5.5 0.825; 5.5 0.875; 5.5 0.925; ...
%     6 0.575; 6 0.86; 6 0.92; 6.5 0.075; 6.5 0.325; 6.5 0.425; 6.5 0.475; 6.5 0.525; 6.5 0.575; 6.5 0.625; 6.5 0.675; 6.5 0.725; 6.5 0.775; ...
%     6.5 0.875; 7 0.575; 7 0.625; 7 0.86; 7.475 0.375; 7.5 0.86; 9.5 0.075; 9.5 0.325; 11.5 0.175;];

fignum = 113;
trouble_region.AR2 = [0.8 Inf];  %if between these rectangular limits, we are in trouble region
factor.trouble = 0.05;  %how much we're allowed to decrease or increase next parameter choice by compared to previous best, if we are in trouble region and new guess is beyond what we've done so far for this body
%also, factor is used to expand lower and upper bounds beyond existing min max over the parameter space for next guess 
factor.regular = 0.2;  %how big/small can amp, lambda, and nlambda grow/shrink compared to most extreme existing values
% factor.bounds = 1.5;  %bounds on optimization guesses, prolly defunct if factor.regular is used
%factor = 1.2 is too big


first_amp = 0.05;  %small amp that should work even for donut bodies


opt_order = {'amp','lambda','nlambda'}; %optimize consecutively in this order

% [AR1_AR2,inds] = unique([ [Results.AR1]' [Results.AR2]' ],'rows');
% AR1_AR2 = setdiff(pms_AR, skips,'rows');



[tails, ~] = unique([ [Results.amp]' [Results.lambda]' [Results.nlambda]'], 'rows');

dep_var = 'Power_eff';
%dep_var  = 'speeds';

max_dep = max([Results.(dep_var)]);
same_cutoff = 1;  %if % diff between all new tail params and best so far is less than this, skip during next runs
clear next_iter

code = [];  param_code = [];

loop = 1:size(AR1_AR2,1);  %unique bodies
%  loop = randperm(size(AR1_AR2,1));

% skip = false(length(loop),1);
new_diffs = NaN(length(loop),3);
%   loop = 252;
%  loop = 286;
% loop = 297;
% loop = 298;

for body_i = loop
    body_i/length(loop)
    body = AR1_AR2(body_i,:);
    
%     if ~isequal(body,[8.5 0.5])   %[7 0.75])
%         continue
%     end

%        0.70343       5.0305       1.7405    0.0059492
%       0.71174       5.1489       1.7973    0.0058559
%       0.82753       5.3745      0.85424    0.0051574
    %  0.48  3.7   1.5
    
    next_iter(body_i).AR1 = body(1);
    next_iter(body_i).AR2 = body(2);
    
    inds = find([Results.AR1] == body(1) & [Results.AR2] == body(2));  %all combos with this body
    
    AR1 = [Results(inds).AR1];
    AR2 = [Results(inds).AR2];
    amp = [Results(inds).amp];  %all amps for this body
    lambda = [Results(inds).lambda]; %all lambdas for this body
    nlambda = [Results(inds).nlambda];  %all nlambdas for this body
    dep = [Results(inds).(dep_var)];
    
    [~, max_ind] = max([Results(inds).(dep_var)]);  %best tail choice for this body
    
    lambda_best = lambda(max_ind);  %best lambda for this body
    nlambda_best = nlambda(max_ind); %best nlambda for this body
    amp_best = amp(max_ind); %best amp for this body
    
    
    amps = amp( lambda == lambda_best & nlambda == nlambda_best);  %all amps that went with best lambda and nlambda
    deps = dep( lambda == lambda_best & nlambda == nlambda_best);  %all effs that went with best lambda and nlambda
    [amps,inds] = sort(amps);  %sorted amps for this body, best lambda and nlambda
    deps = deps(inds);  %now everything is sorted as if plotting eff vs amp
    plotdata.amp = amps;
    plotdata.amp_deps = deps / max_dep;
    
    lambdas = lambda( amp == amp_best & nlambda == nlambda_best);  %all lambdas that went with best amp and nlambda
    deps = dep( amp == amp_best & nlambda == nlambda_best);  %all effs that went with best lambda and nlambda
    [lambdas,inds] = sort(lambdas);  %sorted lambdas for this body, best amp and nlambda
    deps = deps(inds);  %now everything is sorted as if plotting eff vs lambda
    plotdata.lambda = lambdas;
    plotdata.lambda_deps = deps / max_dep;
    
    nlambdas = nlambda( amp == amp_best & lambda == lambda_best);  %all nlambdas that went with best amp and lambda
    deps = dep( amp == amp_best & lambda == lambda_best);
    [nlambdas,inds] = sort(nlambdas);
    deps = deps(inds);  %now everything is sorted as if plotting eff vs amp
    plotdata.nlambda = nlambdas;
    plotdata.nlambda_deps = deps / max_dep;

%     lb = 1/factor*[ min([Results.amp]) min([Results.lambda]) min([Results.nlambda])]';
%     ub = factor*[max([Results.amp]) max([Results.lambda]) max([Results.nlambda])]';
    
        lb = (1-factor.regular)*[ min([amp]) min([lambda]) min([nlambda])]';
    ub = (1+factor.regular)*[max([amp]) max([lambda]) max([nlambda])]';
    
    if isempty(lb) || isempty(ub)  %prolly the first tail for this body
        lb = [ min([Results.amp]) min([Results.lambda]) min([Results.nlambda])]';  %hard limit at tails that worked before - can go beyond this on next sweep
        ub = [max([Results.amp]) max([Results.lambda]) max([Results.nlambda])]';
    end
    
    %       obj1 = @(x) -net([body(1); body(2); x]);
    %   net_guess = fmincon(obj1, [amp_best lambda_best nlambda_best]',[],[],[],[],lb,ub)
    
   pm_ind = find(ismember(pms_AR,body,'rows'));
    pm = pms(pm_ind);
     obj2 = @(x)  - polyvaln(pm,[body(1) body(2) x']);
     
     if isempty(amp_best)  %brand new body - need initial guess....
         amp_best = 0.5;  lambda_best = 3;  nlambda_best = 1.5;
        
     end
     
     pm_guess = fmincon(obj2, [amp_best lambda_best nlambda_best]',[],[],[],[],lb,ub);
     
%      if any(pm_guess < lb) || any(pm_guess > ub)
%          stopafrao
%      end
     
     
     if body(2) >= trouble_region.AR2(1) && body(2) <= trouble_region.AR2(2)  % in trouble region, maybe limit new guess
        %only worried about amp being too big for now
        
        if isempty(amp)  %brand new body in trouble region, hard limit first tail tried
        pm_guess(1) = first_amp;
     
        
        elseif pm_guess(1) > max(amp)  %guessed amp is larger than largest amp tested for this body thus far      don't worry about smaller amps, should be OK
             if pm_guess(1) < max(amp) * (1+factor.trouble)
%                  stopafra
             end
             pm_guess(1) = max(amp) * (1+factor.trouble);
         end
         
     end
         
         next_iter(body_i).param = 'all';
         %     next_iter(body_i).nlambda = net_guess(3);
         %     next_iter(body_i).amp = net_guess(1);
%     next_iter(body_i).lambda = net_guess(2);
    
        next_iter(body_i).nlambda = pm_guess(3);
    next_iter(body_i).amp = pm_guess(1);
    next_iter(body_i).lambda = pm_guess(2);
    
    code(body_i) = NaN;
    param_code(body_i) = NaN;
    
    new_diffs(body_i,:) = [abs(amp_best - pm_guess(1))/amp_best*100    abs(lambda_best - pm_guess(2))/lambda_best*100    abs(nlambda_best - pm_guess(3))/nlambda_best*100 ];
%     if abs(amp_best - pm_guess(1))/amp_best*100 <= same_cutoff && abs(lambda_best - pm_guess(2))/lambda_best*100 <= same_cutoff && abs(nlambda_best - pm_guess(3))/nlambda_best*100 <= same_cutoff
%     skip(body_i) = true;
%     end
    
    % continue   % skip any plotting, just shut yer mouth and run the numbers
      
[amp' lambda' nlambda' dep']
if ~isempty(amp)
pm_dep = polyvaln(pm, [ repmat(body,length(amp),1)  amp' lambda' nlambda']);
else
    pm_dep = [];
end
reldiffs = abs(pm_dep - dep') ./ pm_dep ;
sort(reldiffs)

outlier_cutoff = 0.5;
if ~isempty(dep)
[colors,colorinds] = colordata(200,'jet',[min(dep) max(dep)],dep);
else
    colors = [];  colorinds = [];
end

markersize = 8;
if ~isempty(dep)
figure(830)
% subplot(1,3,1)
for i = 2:length(amp)
plot3(amp(i),lambda(i),nlambda(i),'o','markersize',markersize,'markerfacecolor',colors(colorinds(i),:),'markeredgecolor',colors(colorinds(i),:))
text(amp(i)+( max(amp)-min(amp))*0.02,lambda(i)+( max(lambda)-min(lambda))*0.02,nlambda(i),{num2str(dep(i)/max(dep)),num2str(reldiffs(i))} );
hold on
end
plot3(amp(1),lambda(1),nlambda(1),'p','markersize',markersize+5,'markerfacecolor',colors(colorinds(1),:),'markeredgecolor',colors(colorinds(1),:))
text(amp(1)+( max(amp)-min(amp))*0.02,lambda(1)+( max(lambda)-min(lambda))*0.02,nlambda(1),{num2str(dep(1)/max(dep)),num2str(reldiffs(1))} );
plot3(pm_guess(1),pm_guess(2),pm_guess(3),'*','markersize',markersize+5,'markerfacecolor','k','markeredgecolor','k')
hold off
xlabel('amp'); ylabel('lambda'); zlabel('nlambda');
view([0 90]);
caxis([min(dep) max(dep)])
colormap(colors)
title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))]});
end
% subplot(1,3,2)
% for i = 2:length(amp)
% plot3(amp(i),lambda(i),nlambda(i),'o','markersize',markersize,'markerfacecolor',colors(colorinds(i),:),'markeredgecolor',colors(colorinds(i),:))
% hold on
% end
% plot3(amp(1),lambda(1),nlambda(1),'p','markersize',markersize+5,'markerfacecolor',colors(colorinds(1),:),'markeredgecolor',colors(colorinds(1),:))
% plot3(pm_guess(1),pm_guess(2),pm_guess(3),'*','markersize',markersize+5,'markerfacecolor','k','markeredgecolor','k')
% hold off
% xlabel('amp'); ylabel('lambda'); zlabel('nlambda');
% view([0 0]);
% caxis([min(dep) max(dep)])
% colormap(colors)
% 
% subplot(1,3,3)
% for i = 2:length(amp)
% plot3(amp(i),lambda(i),nlambda(i),'o','markersize',markersize,'markerfacecolor',colors(colorinds(i),:),'markeredgecolor',colors(colorinds(i),:))
% hold on
% end
% plot3(amp(1),lambda(1),nlambda(1),'p','markersize',markersize+5,'markerfacecolor',colors(colorinds(1),:),'markeredgecolor',colors(colorinds(1),:))
% plot3(pm_guess(1),pm_guess(2),pm_guess(3),'*','markersize',markersize+5,'markerfacecolor','k','markeredgecolor','k')
% hold off
% xlabel('amp'); ylabel('lambda'); zlabel('nlambda');
% view([90 0]);
% caxis([min(dep) max(dep)])
% colormap(colors)

% colorbar
      
    %check trend in amp
    diffs = diff(plotdata.amp_deps);
    pause
    continue
    if isempty(diffs) %we only tried one amp so far, try bigger next since curved rods tend to like bigger tails....
        next_iter(body_i).param = 'amp';
        next_iter(body_i).amp = amp_best * factor;
        next_iter(body_i).lambda = lambda_best;
        next_iter(body_i).nlambda = nlambda_best;
        % 'only have one'
        code(body_i) = 0;
        
        figure(fignum)
        
        subplot(1,3,1)
        plot(plotdata.amp,plotdata.amp_deps,'o-');
        ylabel('efficiency / global max efficiency')
        xlabel('tail amplitude (\mum)')
        title(num2str(pm_guess(1)));
        subplot(1,3,2)
        plot(plotdata.lambda,plotdata.lambda_deps,'o-');
        xlabel('tail wavelength (\mum)')
        title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))],num2str(pm_guess(2))});
        subplot(1,3,3)
        plot(plotdata.nlambda,plotdata.nlambda_deps,'o-');
        xlabel('tail # wavelengths')
        title({'only one amp',num2str(pm_guess(3))});
        
        plot_tail_predictions;
        figure(830)
        pause
        
        param_code(body_i) = 1;
    elseif all(diffs > 0)  %eff always increased with amp, try a bigger amp next
        next_iter(body_i).param = 'amp';
        next_iter(body_i).amp = amp_best * factor;
        next_iter(body_i).lambda = lambda_best;
        next_iter(body_i).nlambda = nlambda_best;
        % 'try bigger'
        code(body_i) = 1;
        param_code(body_i) = 1;
        
        figure(fignum)
        
        subplot(1,3,1)
        plot(plotdata.amp,plotdata.amp_deps,'o-');
        ylabel('efficiency / global max efficiency')
        xlabel('tail amplitude (\mum)')
        title(num2str(pm_guess(1)))
        subplot(1,3,2)
        plot(plotdata.lambda,plotdata.lambda_deps,'o-');
        xlabel('tail wavelength (\mum)')
        title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))],num2str(pm_guess(2))});
        subplot(1,3,3)
        plot(plotdata.nlambda,plotdata.nlambda_deps,'o-');
        xlabel('tail # wavelengths')
        title({'try larger amp',num2str(pm_guess(3))})
        
        plot_tail_predictions;
        pause
    elseif all(diffs < 0)  %eff always decreased with amp, try a smaller amp next
        next_iter(body_i).param = 'amp';
        next_iter(body_i).amp = amp_best / factor;
        next_iter(body_i).lambda = lambda_best;
        next_iter(body_i).nlambda = nlambda_best;
        %  'try smaller'
        code(body_i) = -1;
        param_code(body_i) = 1;
        
        figure(fignum)
        
        subplot(1,3,1)
        plot(plotdata.amp,plotdata.amp_deps,'o-');
        ylabel('efficiency / global max efficiency')
        xlabel('tail amplitude (\mum)')
        title(num2str(pm_guess(1)))
        subplot(1,3,2)
        plot(plotdata.lambda,plotdata.lambda_deps,'o-');
        xlabel('tail wavelength (\mum)')
        title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))],num2str(pm_guess(2))});
        subplot(1,3,3)
        plot(plotdata.nlambda,plotdata.nlambda_deps,'o-');
        xlabel('tail # wavelengths')
        title({'try smaller amp',num2str(pm_guess(3))})
        
        plot_tail_predictions;
        pause
        
    else %there must be at least one rel min or max
        sign_changes = diff(sign(diffs));  %basically concavity, should be zeros for monotonic segments, and a +2 for a each rel min, a -2 for each single rel max
        sign_changes(sign_changes == 0) = []; %get rid of monotonic sequences, leaving any rel mins, maxes
        
        if sign_changes == 2
            disp('found relative min for amp')
            pause
        elseif sign_changes == -2  %there is one rel max, we're done for amp
              
            %next, look at lambda
                 
            diffs = diff(plotdata.lambda_deps);
                      
            if isempty(diffs) %we only tried one lambda so far, try bigger next since curved rods tend to like bigger tails....
                next_iter(body_i).param = 'lambda';
                next_iter(body_i).lambda = lambda_best * factor;
                next_iter(body_i).amp = amp_best;
                next_iter(body_i).nlambda = nlambda_best;
                % 'only have one'
                code(body_i) = 0;
                param_code(body_i) = 2;
                
                figure(fignum)
                
                subplot(1,3,1)
                plot(plotdata.amp,plotdata.amp_deps,'o-');
                ylabel('efficiency / global max efficiency')
                xlabel('tail amplitude (\mum)')
                title(num2str(pm_guess(1)))
                subplot(1,3,2)
                plot(plotdata.lambda,plotdata.lambda_deps,'o-');
                xlabel('tail wavelength (\mum)')
                title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))],num2str(pm_guess(2))});
                subplot(1,3,3)
                plot(plotdata.nlambda,plotdata.nlambda_deps,'o-');
                xlabel('tail # wavelengths')
                title({'only one lambda',num2str(pm_guess(3))})
                
                plot_tail_predictions;
                pause
                
            elseif all(diffs > 0)  %eff always increased with lambda, try a bigger lambda next
                next_iter(body_i).param = 'lambda';
                next_iter(body_i).lambda = lambda_best * factor;
                next_iter(body_i).amp = amp_best;
                next_iter(body_i).nlambda = nlambda_best;
                % 'try bigger'
                code(body_i) = 1;
                param_code(body_i) = 2;
                
                figure(fignum)
                
                subplot(1,3,1)
                plot(plotdata.amp,plotdata.amp_deps,'o-');
                ylabel('efficiency / global max efficiency')
                xlabel('tail amplitude (\mum)')
                title(num2str(pm_guess(1)))
                subplot(1,3,2)
                plot(plotdata.lambda,plotdata.lambda_deps,'o-');
                xlabel('tail wavelength (\mum)')
                title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))],num2str(pm_guess(2))});
                subplot(1,3,3)
                plot(plotdata.nlambda,plotdata.nlambda_deps,'o-');
                xlabel('tail # wavelengths')
                title({'try larger lambda',num2str(pm_guess(3))})
                
                plot_tail_predictions;
                pause
            elseif all(diffs < 0)  %eff always decreased with lambda, try a smaller lambda next
                next_iter(body_i).param = 'lambda';
                next_iter(body_i).lambda = lambda_best / factor;
                next_iter(body_i).amp = amp_best;
                next_iter(body_i).nlambda = nlambda_best;
                %  'try smaller'
                code(body_i) = -1;
                param_code(body_i) = 2;
                %pause
                
                figure(fignum)
                
                subplot(1,3,1)
                plot(plotdata.amp,plotdata.amp_deps,'o-');
                ylabel('efficiency / global max efficiency')
                xlabel('tail amplitude (\mum)')
                title(num2str(pm_guess(1)))
                subplot(1,3,2)
                plot(plotdata.lambda,plotdata.lambda_deps,'o-');
                xlabel('tail wavelength (\mum)')
                title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))],num2str(pm_guess(2))});
                subplot(1,3,3)
                plot(plotdata.nlambda,plotdata.nlambda_deps,'o-');
                xlabel('tail # wavelengths')
                title({'try smaller lambda',num2str(pm_guess(3))})
                
                plot_tail_predictions;
                pause
            else %there must be a rel min or max for lambda
                sign_changes = diff(sign(diffs));  %basically concavity, should be zeros for monotonic segments, and a +2 for a each rel min, a -2 for each single rel max
                sign_changes(sign_changes == 0) = []; %get rid of monotonic sequences, leaving any rel mins, maxes
                
                if sign_changes == 2
                    disp('found relative min for lambda')
                    pause
                elseif sign_changes == -2  %there is one rel max, we're done for lambda
                    %             next_iter(body_i).param = NaN;
                    %             next_iter(body_i).lambda = NaN;
                    %                     next_iter(body_i).amp = NaN;
                    %         next_iter(body_i).nlambda = NaN;
                    %             code(body_i) = NaN;
                    
                    
                    
                    diffs = diff(plotdata.nlambda_deps);
                    
                    
                    if isempty(diffs) %we only tried one nlambda so far, try bigger next since curved rods tend to like bigger tails....
                        next_iter(body_i).param = 'nlambda';
                        next_iter(body_i).nlambda = nlambda_best * factor;
                        next_iter(body_i).amp = amp_best;
                        next_iter(body_i).lambda = lambda_best;
                        % 'only have one'
                        code(body_i) = 0;
                        param_code(body_i) = 3;
                        
                        figure(fignum)
                        
                        subplot(1,3,1)
                        plot(plotdata.amp,plotdata.amp_deps,'o-');
                        ylabel('efficiency / global max efficiency')
                        xlabel('tail amplitude (\mum)')
                        title(num2str(pm_guess(1)))
                        subplot(1,3,2)
                        plot(plotdata.lambda,plotdata.lambda_deps,'o-');
                        xlabel('tail wavelength (\mum)')
                        title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))],num2str(pm_guess(2))});
                        subplot(1,3,3)
                        plot(plotdata.nlambda,plotdata.nlambda_deps,'o-');
                        xlabel('tail # wavelengths')
                        title({'only one nlambda',num2str(pm_guess(3))})
                        plot_tail_predictions;
                        pause
                    elseif all(diffs > 0)  %eff always increased with nlambda, try a bigger nlambda next
                        next_iter(body_i).param = 'nlambda';
                        next_iter(body_i).nlambda = nlambda_best * factor;
                        next_iter(body_i).amp = amp_best;
                        next_iter(body_i).lambda = lambda_best;
                        % 'try bigger'
                        code(body_i) = 1;
                        param_code(body_i) = 3;
                        
                        figure(fignum)
                        
                        subplot(1,3,1)
                        plot(plotdata.amp,plotdata.amp_deps,'o-');
                        ylabel('efficiency / global max efficiency')
                        xlabel('tail amplitude (\mum)')
                        title(num2str(pm_guess(1)))
                        subplot(1,3,2)
                        plot(plotdata.lambda,plotdata.lambda_deps,'o-');
                        xlabel('tail wavelength (\mum)')
                        title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))],num2str(pm_guess(2))});
                        subplot(1,3,3)
                        plot(plotdata.nlambda,plotdata.nlambda_deps,'o-');
                        xlabel('tail # wavelengths')
                        title({'try larger nlambda',num2str(pm_guess(3))})
                        plot_tail_predictions;
                        pause
                    elseif all(diffs < 0)  %eff always decreased with nlambda, try a smaller nlambda next
                        next_iter(body_i).param = 'nlambda';
                        next_iter(body_i).nlambda = nlambda_best / factor;
                        next_iter(body_i).amp = amp_best;
                        next_iter(body_i).lambda = lambda_best;
                        %  'try smaller'
                        code(body_i) = -1;
                        param_code(body_i) = 3;
                        %pause
                        
                        figure(fignum)
                        
                        subplot(1,3,1)
                        plot(plotdata.amp,plotdata.amp_deps,'o-');
                        ylabel('efficiency / global max efficiency')
                        xlabel('tail amplitude (\mum)')
                        title(num2str(pm_guess(1)))
                        subplot(1,3,2)
                        plot(plotdata.lambda,plotdata.lambda_deps,'o-');
                        xlabel('tail wavelength (\mum)')
                        title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))],num2str(pm_guess(2))});
                        subplot(1,3,3)
                        plot(plotdata.nlambda,plotdata.nlambda_deps,'o-');
                        xlabel('tail # wavelengths')
                        title({'try smaller nlambda',num2str(pm_guess(3))})
                        plot_tail_predictions;
                        pause
                        
                    else %there must be a rel min or max
                        sign_changes = diff(sign(diffs));  %basically concavity, should be zeros for monotonic segments, and a +2 for a each rel min, a -2 for each single rel max
                        sign_changes(sign_changes == 0) = []; %get rid of monotonic sequences, leaving any rel mins, maxes
                        
                        if sign_changes == 2
                            disp('found relative min for nlambda')
                            pause
                        elseif sign_changes == -2  %there is one rel max, we're done for nlambda
                            next_iter(body_i).param = NaN;
                            next_iter(body_i).lambda = NaN;
                            next_iter(body_i).amp = NaN;
                            next_iter(body_i).nlambda = NaN;
                            code(body_i) = NaN;
                            param_code(body_i) = NaN;
                            %
                            figure(fignum)
                            
                            subplot(1,3,1)
                            plot(plotdata.amp,plotdata.amp_deps,'o-');
                            %                             hold on
                            %                                                         coeffs = polyfit(plotdata.amp,plotdata.amp_deps, 2);
                            %                             fit_x = linspace(min(plotdata.amp),max(plotdata.amp), 100);
                            %                             fit_y = polyval(coeffs,fit_x);
                            %                             plot(fit_x,fit_y, '-','linewidth',2)
                            %                             hold off
                            ylabel('efficiency / global max efficiency')
                            xlabel('tail amplitude (\mum)')
                            title(num2str(pm_guess(1)))
                            
                            subplot(1,3,2)
                            plot(plotdata.lambda,plotdata.lambda_deps,'o-');
                            %                                       hold on
                            %                                                         coeffs = polyfit(plotdata.lambda,plotdata.lambda_deps, 2);
                            %                             fit_x = linspace(min(plotdata.lambda),max(plotdata.lambda), 100);
                            %                             fit_y = polyval(coeffs,fit_x);
                            %                             plot(fit_x,fit_y, '-','linewidth',2)
                            %                             hold off
                            xlabel('tail wavelength (\mum)')
                            title({['AR1 ',num2str(body(1)),'   AR2 ',num2str(body(2))],num2str(pm_guess(2))});
                            
                            
                            subplot(1,3,3)
                            plot(plotdata.nlambda,plotdata.nlambda_deps,'o-');
                            %                                       hold on
                            %                                                         coeffs = polyfit(plotdata.nlambda,plotdata.nlambda_deps, 2);
                            %                             fit_x = linspace(min(plotdata.nlambda),max(plotdata.nlambda), 100);
                            %                             fit_y = polyval(coeffs,fit_x);
                            %                             plot(fit_x,fit_y, '-','linewidth',2)
                            %                             hold off
                            
                            xlabel('tail # wavelengths')
                            title({'Everything optimized',num2str(pm_guess(3))})
                            
                            plot_tail_predictions;
                            drawnow
                            
                            pause
                            
                            
                            
                            
                            
                        else %there is more than one rel min or max for nlambda....
                            disp('more than one rel min or max for nlambda!')
                            pause
                        end % various rel min, max situations for nlambda
                        
                    end %various possiblities for nlambda
                    
                else %there is more than one rel min or max for lambda....
                    figure(456)
                    plot(lambdas,deps,'o-'); grid on;
                    disp('more than one rel min or max for lambda!')
                    pause
                end
                
            end
            
        else %there is more than one rel min or max for amp....
            disp('more than one rel min or max for amp!')
            pause
        end
        
    end
    
    
end


temp = [next_iter.amp];  inds = ~isnan(temp);  %get rid of combos that are already halfway decently optimized
new_sweep.AR1 = [next_iter(inds).AR1];
new_sweep.AR2 = [next_iter(inds).AR2];
new_sweep.amp = [next_iter(inds).amp];
new_sweep.lambda = [next_iter(inds).lambda];
new_sweep.nlambda = [next_iter(inds).nlambda];

optimized_frac = sum(isnan([next_iter.amp])) / numel([next_iter.amp])

new_sweep_orig = new_sweep;
sstopoo
%%
new_sweep = new_sweep_orig;
tol = 0;  %don't bother with runs with predicted improvement less than tol %
skip = all(new_diffs < tol, 2);

new_sweep.AR1 = new_sweep.AR1(~skip);
new_sweep.AR2 = new_sweep.AR2(~skip);
new_sweep.amp = new_sweep.amp(~skip);
new_sweep.lambda = new_sweep.lambda(~skip);
new_sweep.nlambda = new_sweep.nlambda(~skip);

randomize_new_sweep;  % randomizes order of new_sweep runs in-place

new_sweep

%%
V = 1;  sphererad = (V*3/4/pi)^(1/3);  %equivalent sphere radius

new_tails = unique([  [next_iter.amp];  [next_iter.lambda];  [next_iter.nlambda]  ]','rows');
new_tails = new_tails( ~all(isnan(new_tails),2),:) ;  %get rid of new_iters that actually are already OK for all 3 params

%             amp        lambda      nlambda
base_tails = [0.402  4.68*sphererad  1.49];

new_tails ./ repmat(base_tails,size(new_tails,1),1);

%new_tails = setdiff(new_tails, tails,'rows') ; %don't bother recreating tail meshes we already have

randomize_new_sweep;

save('E:\Hull\sweeps\new_sweep.mat','new_sweep');

stopa
new_sweep2 = new_sweep;
load('E:\Hull\sweeps\new_sweep.mat','new_sweep');

shat = setdiff([new_sweep2.AR1' new_sweep2.AR2' new_sweep2.amp' new_sweep2.lambda' new_sweep2.nlambda'],[new_sweep.AR1' new_sweep.AR2' new_sweep.amp' new_sweep.lambda' new_sweep.nlambda'], 'rows')

new_sweep.AR1 = shat(:,1)'; new_sweep.AR2 = shat(:,2)';
new_sweep.amp = shat(:,3)'; new_sweep.lambda = shat(:,4)'; new_sweep.nlambda = shat(:,5)';

randomize_new_sweep;
%repmat(base_tails,size(new_tails,1),1) ./  new_tails