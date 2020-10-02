clf(fignum) 

labels = {  '{\itU}'  ,  '{\it\Omega}^B'  ,  '{\it\Omega}^F'};
fontsize = 16;

c = 0; current_diffs = [];  pc = 0; ylims = [];
    for f = 1:length(interp_fields)
        if input.accuracy.interpolation.vector_normalization % usual vector magnitude for U and Omega.  in case of bacteria, vector magnitude of tail rotation rate omega is still omega, so no change there
            vector_mag = sqrt( sum( Solutions.(interp_fields{f}).^2 , 2 ) ); %vector magnitude at each phase
            interp_mag = trig_interp_fit(Solutions.phase, vector_mag);
        end
        
        for i = 1:size(Solutions.(interp_fields{f}),2) %scalar components of this variable, e.g. U(1:3), Omega(1:3).  trig_interp only works for scalars currently
            
            interpolant(n_iter).(interp_fields{f})(i) = trig_interp_fit(Solutions.phase,Solutions.(interp_fields{f})(:,i));
            
            if n_iter >= 2 %we have enough iterations to compare the last two
                fun2 = @(x) trig_interp_eval(interpolant(n_iter).(interp_fields{f})(i), x);   %more accurate
                fun1 = @(x) trig_interp_eval(interpolant(n_iter-1).(interp_fields{f})(i), x); %less accurate
                if input.accuracy.interpolation.vector_normalization
                    fun_norm = @(x) trig_interp_eval(interp_mag,x);
                else
                    fun_norm = fun_2;
                end
                
                rel_diff(n_iter).(interp_fields{f})(i) = RMS_fun_diff(fun1, fun2, fun_norm, phase_bounds);  %global phase angle should always be defined from 0 - 2*pi, and will contain entire cyclic motion
                c = c+1;
                current_diffs(c) = rel_diff(n_iter).(interp_fields{f})(i);  %concatenate all rel_diffs in a vector for convenience below
            else
                rel_diff(n_iter).(interp_fields{f})(i) = NaN;  %placeholder for first undefined relative difference
            end
            
            if input.output.interpolation.doplot  %plot interpolant curves
                figure(fignum)
                pc = pc+1;
                subplot(numplots,1,pc)
                cla
            
                
                  if n_iter >= 3 %plot 2nd to last solutions if they exist
                    plot(x, trig_interp_eval(interpolant(n_iter-2).(interp_fields{f})(i), x),'b--','linewidth',1.5); hold on;
%                     plot([Solutions.phase(1:4:end); phase_bounds(2)],[Solutions.(interp_fields{f})(1:4:end,i); Solutions.(interp_fields{f})(1,i)] ,'bo','markersize',18);
                     plot([Solutions.phase(1:4:end);],[Solutions.(interp_fields{f})(1:4:end,i); ] ,'bo','markersize',18);
                  end
                  
                    if n_iter >= 2 %plot 2nd to last solutions if they exist
                    plot(x, trig_interp_eval(interpolant(n_iter-1).(interp_fields{f})(i), x),'r--','linewidth',1.5); hold on;
%                     plot([Solutions.phase(1:2:end); phase_bounds(2)],[Solutions.(interp_fields{f})(1:2:end,i); Solutions.(interp_fields{f})(1,i)] ,'ro','markersize',10);
                     plot([Solutions.phase(1:2:end); ],[Solutions.(interp_fields{f})(1:2:end,i);] ,'ro','markersize',10);
                end
                  
                plot(x, trig_interp_eval(interpolant(n_iter).(interp_fields{f})(i), x),'k-','linewidth',1.5); hold on;
%                 plot([Solutions.phase(:); phase_bounds(2)],[Solutions.(interp_fields{f})(:,i); Solutions.(interp_fields{f})(1,i)],'ko','markersize',4,'markerfacecolor','k');
                 plot([Solutions.phase(:);],[Solutions.(interp_fields{f})(:,i);],'ko','markersize',4,'markerfacecolor','k');
                xlim(phase_bounds)
                ylim auto
                grid on
                set(gca,'fontsize',13);
                if n_iter >= 2
%                     title(['relative difference = ',num2str(current_diffs(c))]);
                end
                if size(Solutions.(interp_fields{f}),2) > 1
                    ylabel([labels{f},'_',num2str(i)],'fontsize',fontsize);
                else
                    ylabel([labels{f}],'fontsize',fontsize);
                end
                ylims{f}(i,:) = ylim;
                
            end
            
            
        end
    end
    
    xlabel( ['{\it\theta}^{',char(8201),'F} (rad)'],'fontsize',fontsize);