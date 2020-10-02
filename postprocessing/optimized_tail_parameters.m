%% Plot interpolants for optimized amp, lambda, nlambda for Dm, eff, temporal

cvals = { [],[],[],[],[],[]};
clear text
%
%  [0.85 0.9 0.95 0.97 0.98 0.99 1  1.006 1.0075 1.008 1.009 ],... adj speed

npts = 150;

fullnames = {};  metrics = {}; displaynames = {};
% try, delete(handles); end
Fs_tails = {};  handles = [];
% options = optimset('MaxFunEvals',1E9,'MaxIter',1E9,'TolX',1E-6,'TolFun',1E-8,'display','iter');
options = optimoptions('patternsearch','display','off','FunctionTolerance',1E-8,'MaxFunctionEvaluations',1E9,'MaxIterations',1E9,'MeshTolerance',1E-8,'StepTolerance',1E-8);
normalize_metrics = false; % normalize to metric value for a sphere

% fields = {'Dm','eff','temporal'};  display_fields = {'Dispersal','Swimming Efficiency','Temporal S/N'}; file_fields = {'Dispersal','Efficiency','Temporal_SN'};
fields = {'eff','temporal'};   % fields = {'eff'};

display_fields = {'\Psi_{swim}','\Psi_{chemo}'}; file_fields = {'Efficiency','Temporal_SN'};
parameters = {'amp','lambda','nlambda'};  display_parameters = {'\ita\rm','\lambda','\itn_{\lambda}'};  file_parameters = {'amplitude','wavelength','num_wavelengths'};

panel_labels = {'A','B','C','D','E','F'};
clear all_Z
fignum = 0;
sp = 0;
for p = 1:length(parameters)
    
    
    limits = [1 10; 0 1];
    
    [X,Y] = meshgrid(linspace(limits(1,1),limits(1,2),npts),linspace(limits(2,1),limits(2,2),npts));
    [standardized_XY] = standardize([X(:) Y(:)], limits);
    standardized_X = NaN(size(X));  standardized_X(:) = standardized_XY(:,1);
    standardized_Y = NaN(size(Y));  standardized_Y(:) = standardized_XY(:,2);
    
    fontsizes.labels = 16;
    
    for m = 1:length(fields)
        sp = sp+1;
        if ~normalize_metrics
            Z = F.(fields{m}).F.(parameters{p})(standardized_X,standardized_Y);
        else
            Z = F.(fields{m}).F.(parameters{p})(standardized_X,standardized_Y) / F.(fields{m}).F.(parameters{p})(0,0);
        end
        
        if p == 3 && m == 1
            Z = all_Z(:,:,1) ./ all_Z(:,:,3);
        end
        
        figure(328); set(gcf,'Color','w');
%         subtightplot(3,2,sp,[0.02 0.000]);
        subaxis(3,2,sp,'spacing',0.03)
        
        %         fignum = fignum + 1;
        %         handles(m) = figure(fignum);  set(handles(m),'Name',[ fields{m}  , '  ', parameters{p} ]);
        
        all_Z(:,:,sp) = Z;
        
        pcolor(X,Y,Z);  shading interp;  cb = colorbar;
        
        units = {'(\mum)' , '(\mum)' , ''};
        if ismember(sp, [2 4 6])
          cb.Label.String = [display_parameters{p},' ',units{p}];
        end
          cb.FontSize = 12;
        cb.Label.FontSize = 16;
      
        
        %                 title([ display_fields{m}  , ' optimized ', display_parameters{p} ],'interpreter','none');
                if ismember(sp,[1 2])
                    tit =  title([ display_fields{m}  , ' - optimized flagella']); tit.FontSize = 12;
                end
                
                
%                  title([ display_parameters{p} ],'interpreter','none');
          
        
%         grid on
        hold on
        
%         if ~isempty(cvals{m})
%             [C,ch] = contour(unique(X),unique(Y),Z,cvals{m},'k--','linewidth',1);
%         else
%             [C,ch] = contour(unique(X),unique(Y),Z,20,'k--','linewidth',1);
%         end
%         clabel(C,ch,'fontsize',12,'LabelSpacing',1000);
        
        
        
        if isfield(AR1_AR2_data,fields{m})
             d = plot(AR1_AR2_data.(fields{m})(:,1),AR1_AR2_data.(fields{m})(:,2), '.','markersize',2,'markerfacecolor','none','markeredgecolor',repmat(0.4,1,3));
        end
        
        fun = @(x) - F.(fields{m}).F.(parameters{p})(x(1),x(2));
        guesses = [1 0; 4 0.7; 10 0.9;  10 0;  5.5 0.5;  1.5 0];
        [standardized_guesses] = standardize(guesses, limits);
        clear body metric_val
        for g = 1:size(guesses,1)
            %         [body(g,:), metric_val(g)] = fminsearch(fun, standardized_guesses(g,:)',options);
            try
            [body(g,:), metric_val(g)]  = patternsearch(fun, standardized_guesses(g,:)',[],[],[],[],[0 0]',[1 1]',[],options);
            catch
                body(g,:) = [NaN NaN];  metric_val(g) = NaN;
            end
        end
        [~,ind] = max(-metric_val);  best_body = unstandardize(body(ind,:),limits);
        
%         switch fields{m}
%             case {'eff','construction','tumbling'}
%                 starh = plot(best_body(1),best_body(2),'kp','markersize',18,'markerfacecolor','k');
%             otherwise
%                 temp = 0.3;
%                 starh = plot(best_body(1),best_body(2),'p','markersize',18,'markeredgecolor',[temp temp temp],'markerfacecolor',[temp temp temp]);
%         end
        
        
        hold off
        set(gca,'fontsize',12);
        if sp == 5
        xlabel('$$ \mathcal{L} $$ (elongation)','interpreter','latex','fontsize',17)
        ylabel('$$ \mathcal{K} $$ (elongation)','interpreter','latex','fontsize',17)
        end
        
        if ismember(sp,[1 2 3 4 6])
            set(gca,'XTickLabel',[]);  set(gca,'YTickLabel',[]);
        end
        if ismember(sp,[5])
            set(gca,'XTick',[1 2:1:10]);  
        end
  
      set(gca,'PlotBoxAspectRatio',[   1      0.80206      0.80206]);
   
      set(gcf,'Position',[ -1289          42         992         954]);
        %     pause
        %     switch metric
        %         case 'construction'
        %             saveas(handles(m),['C:\Hull\pareto\tasks\',fullnames{m},'.png']);
        %         otherwise
        %             saveas(handles(m),['C:\Hull\pareto\tasks\',displaynames{m},'.png']);
        %     end
%         saveas(gcf, ['C:\Hull\current graphs\optimized tail parameters\', file_fields{m}  , ' optimized ', file_parameters{p},'.fig' ]);
%         saveas(gcf, ['C:\Hull\current graphs\optimized tail parameters\',  file_fields{m}  , ' optimized ', file_parameters{p},'.png' ]);
        text(1.1,0.95,panel_labels{sp},'FontWeight','bold','FontSize',14);
    end
    
end
return
%%
col = @(M) M(:);
clear tail;  tail(:,1) = col(all_Z(:,:,1));  tail(:,2) = col(all_Z(:,:,3));  tail(:,3) = col(all_Z(:,:,5));
[arclength] = bacterial_tail_arclength(tail);

