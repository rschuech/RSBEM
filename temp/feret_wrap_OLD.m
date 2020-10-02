% widths = test_data.width;
% curvatures = test_data.curvature;
% Ferets = test_data.length;
% Lengths = NaN(size(widths));

do_width_correction = true;  % convert measured mean width to geometric curved rod width
do_curv_correction = true;  % convert measured curvature2 based on pole points to correct curvature
% do_width_correction = false;


do_plots = false;
% do_plots = true;


% bad_inds = [];

parfor im = 1:length(data) % images
    
    
    
    if ~isfield(data(im),'length_fixed') || isempty(data(im).length_fixed)
        data(im).length_fixed = NaN(size(data(im).length));
        data(im).width_fixed = NaN(size(data(im).width));
        data(im).curvature_fixed = NaN(size(data(im).curvature));
    end
    
    
    donuts = false(1,length(data(im).length));
    Curv_sol = []; u_sol_best = [];
    for c = 1:length(data(im).length)  %individual cells
        c
        if ~isnan( data(im).length_fixed(c) )
                        continue
        end
        
        
        
        donut_safety_factor = 0.05;  % if known Feret value is within this factor of the Feret of a semicircular-or-greater rod, then switch to PoleDist method
        SF1_max = 50; % no way SF1 can ever be larger than this
        
        constants = [];
        constants.width = data(im).width(c);
        
        
        constants.curvature = data(im).curvature(c);
        
        Feret = data(im).length(c);
        PoleDist = data(im).poledist(c);
        
        PoleDist = NaN;
        
        constants.curvature = max(constants.curvature , 1E-9); % set curvature < 1E-10 to 1E-10 to avoid numerical precision problem whereby the code breaks for very straight rods
        global_min_Length = constants.width;  % if shape is a circle
        global_max_Length = min( 2*pi*(1/constants.curvature + constants.width/2)  ,  SF1_max*constants.width); % circumference of closed circular centerline
        
        %% determine min and max geometrically possible Lengths such that curved rod is non self intersecting
        %  options = optimoptions('fmincon', 'useparallel',false,'display','final');
        options = optimoptions('patternsearch','MeshTolerance',1E-8,'ConstraintTolerance',1E-8,'display','off');
        nonlincon = @(Length) is_curved_rod_self_intersected(Length, constants, do_width_correction);
        
%         objfun = @(Length) Length;
%         %         Length_guess = (global_min_Length + global_max_Length) / 2;
%         Length_guess = global_min_Length * 1.05;
%         % min_Length = fmincon(objfun,Length_guess, [],[],[],[],global_min_Length, global_max_Length, nonlincon,options);
%         [min_Length,~,flag] = patternsearch(objfun,Length_guess, [],[],[],[],global_min_Length, global_max_Length, nonlincon,options);
%         if flag <= 0
%             error('min_Length calc failed');
%         end
        
        min_Length = constants.width;
        min_Feret = Length2Feret(min_Length, constants,true, do_width_correction);
        
        
        
        
        
        objfun = @(Length) -Length;
        %         Length_guess = (global_min_Length + global_max_Length) / 2;
        Length_guess = global_min_Length * 1;
        % max_Length = fmincon(objfun,Length_guess, [],[],[],[],global_min_Length, global_max_Length, nonlincon,options);
        [max_Length,~,flag] = patternsearch(objfun,Length_guess, [],[],[],[],global_min_Length, global_max_Length, nonlincon,options);
        if flag <= 0
            error('max_Length calc failed');
        end
        max_Feret = Length2Feret(max_Length, constants,true, do_width_correction);
        
        %%
        
        
        if Feret < (1 - donut_safety_factor) * max_Feret
            donut = false;
            if ~do_curv_correction
                
                
                Length_guess = (min_Length + max_Length) / 2;
                
                
                objfun = @(Length) ( Feret - Length2Feret(Length, constants, false, do_width_correction) );
                
                % [Length_sol, fval] = fminbnd(objfun , min_Length, max_Length,optimset('display','iter'))
                
                [Length_sol, fval] = fzero(objfun , [min_Length max_Length]);
                if abs(fval) > 1E-6
                    error('fzero failed');
                end
                
                
                
                
            else
                
                %%
                Feret_Curv2 = [Feret , constants.curvature];
                check_intersections = true;  % makes almost no difference in speed so safer to do it
                objfun = @(Length_Curv) objfun_feret_curv(Length_Curv, Feret_Curv2, constants, check_intersections, do_width_correction);
                lb = [min_Length , constants.curvature * 1E-1];  ub = [max_Length , constants.curvature * 3];
                %             guess = [Length_sol, constants.curvature];
                guess = [Feret, constants.curvature* 1];  % * 3
                %                 guess = [258, constants.curvature* 1.2];  % * 3
                
                %                 options = optimoptions('patternsearch','StepTolerance',1E-6/1,'MeshTolerance',1E-6/1,'FunctionTolerance',1E-12,'MaxIterations',1E6,'ConstraintTolerance',1E-8,'display','iter');
                
                %                 [Length_Curv,fval,flag] = patternsearch(objfun,guess, [],[],[],[],lb, ub, [],options);
                
                options = optimoptions('fsolve','FunctionTolerance',1E-10, 'MaxFunctionEvaluations', 1E5, 'MaxIterations',1E4,'display','off');
                [Length_Curv,fval,flag] = fsolve(objfun, guess,options);
                
                
                
                if flag <= 0
                    error('Length and Curvature fix calc completely failed');
                end
                
                [~, Feret_sol, Curv2_sol] = objfun(Length_Curv);
                ansy = [Feret_Curv2;  Feret_sol Curv2_sol];
                
                abs_errors = abs( [Feret - Feret_sol , constants.curvature - Curv2_sol ] );
                
                if any(abs_errors > 1E-6)
                    error('Abs errors too big');
                end
                
                Length_sol = Length_Curv(1);
                Curv_sol = Length_Curv(2);
                
                
                 width = 2*Length_sol/(pi - 4)*( sqrt( 1 + (pi - 4)/Length_sol*constants.width ) - 1 );
                  
                 shaperoo = [  Length_sol / width   ,     Length_sol / (2*pi*1/Curv_sol) ]
                
                % computed minlength 131.93    fval 7.86E-18    length
                % 133.87
                
                 % computed minlength 103.62    fval 3.55E-18    length
                % 127.48
                
            end
            constants2 = constants;  constants2.curvature = Curv_sol;
            [~,u_sol_best] = Length2Feret(Length_sol, constants2, true, do_width_correction);
            
        else
            donut = true;
            donuts(c) = true;
            disp(['Feret cutoff exceeded for ',data(im).file,'    ','cell index ',num2str(c)]);
%             bad_inds = unique([bad_inds   im]);
            
            if isnan(PoleDist)
                Length_sol = NaN;
                Curv_sol = NaN;
                disp(['No PoleDist for ',data(im).file,'    ','cell index ',num2str(c)]);
                continue
            end
            objfun = @(Length) - Length2PoleDist(Length, constants, false, do_width_correction);
            Length_guess = (min_Length + max_Length)/2;
            [semicircular_Length] = patternsearch(objfun,Length_guess, [],[],[],[],min_Length, max_Length, [],options);  %pole distance for approx. semicircular rod
            
            objfun = @(Length) PoleDist - Length2PoleDist(Length, constants, false, do_width_correction);
            [Length_sol2, fval2] = fzero(objfun , [semicircular_Length    max_Length]);  % there are usually 2 possible lengths for the same PoleDist; we want the larger one corresponding to a donut-like rod
            Feret2 = Length2Feret(Length_sol2, constants,true, do_width_correction);
            
            if sign(objfun(min_Length)) == sign(objfun(semicircular_Length))  % the less curved rod solution is impossible
                Length_sol1 = NaN;  fval1 = NaN;  Feret1 = NaN;
            else
                [Length_sol1, fval1] = fzero(objfun , [min_Length   semicircular_Length  ]);  % smaller length solution going with less curved rod
                Feret1 = Length2Feret(Length_sol1, constants,true, do_width_correction);
            end
            
            if abs(fval2) > 1E-6  ||  abs(fval1) > 1E-6
                error('fzero failed');
            end
            
            
            
            
            if abs( (Length_sol1 - Length_sol2) / min([Length_sol1 Length_sol2])) < 0.05  % which curved rod is the correct one is hard to tell, it is right on the tipping point?
                disp('shateroo')
                pause
            end
            
            if abs(Feret - Feret1) < abs(Feret - Feret2)
                Length_sol = Length_sol1;
            else
                Length_sol = Length_sol2;
            end
            
            %              Length_sol = Length_sol1;
            Length_sol = Length_sol2;
            
        end
        
        if do_width_correction
            constants.mean_width = constants.width; % original measured value
            constants.width = 2*Length_sol/(pi - 4)*( sqrt( 1 + (pi - 4)/Length_sol*constants.width ) - 1 );  % from mean width to actual width
        end
        
        if do_curv_correction
            constants.curvature = Curv_sol;
        end
        
        if do_plots
            constants.length = Length_sol - constants.width;
            [~,~,~,~,u1,u2,u3,shift] = curved_rod_parameters(constants);  u0 = 0;
            u = unique( [ linspace(0,1,500)   u0 u1 u2 u3 1 ] );
            [x,y] = curved_rod_pts(u,constants,shift);
            
            figure(45);  plot(x,y,'-','linewidth',3);  hold on
            plot( [ x(u == u0) x(u == u1) ] ,  [ y(u == u0) y(u == u1) ] ,'k:','linewidth',1.5);
            plot( [ x(u == u2) x(u == u3) ] ,  [ y(u == u2) y(u == u3) ] ,'k:','linewidth',1.5);
            if ~donut
                [x,y] = curved_rod_pts(u_sol_best,constants,shift);
                plot(x,y,'k--o','linewidth',1.5,'markerfacecolor','k');
            end
            
            cx = 0 + constants.width/2;
            cy = 0 + 1/constants.curvature;
            theta_range = (Length_sol - constants.width) / (1/constants.curvature);
            theta_min = -pi/2;  theta_max = theta_min + theta_range;
            
            
            t = linspace(theta_min,theta_max,100);
            x_arc = cx + 1/constants.curvature * cos(t);
            y_arc = cy + 1/constants.curvature * sin(t);
            [xpole ,ypole]= curved_rod_pts( (u2+u3)/2,constants,shift); % final pole
            [xtemp1,ytemp1] = curved_rod_pts(u2,constants,shift);
            [xtemp2,ytemp2] = curved_rod_pts(u3,constants,shift);
            xmid = (xtemp1+xtemp2)/2; % center of final circle
            ymid = (ytemp1+ytemp2)/2;
            
            x_centerline = [0 (0 + constants.width/2)  x_arc  xmid  xpole];
            y_centerline = [0            0             y_arc  ymid  ypole];
            
            plot(x_centerline,y_centerline,'r--','linewidth',1.5);
            
            hold off
            axis equal
            grid on
            %     title(['ind = ',num2str(L),'     ','feret = ',num2str(feret_best)]);
            %         title([  'ind = ',num2str(i),'    ','SF1 = ',num2str(actual(i,1)),'    ','SF2 = ',num2str(actual(i,2))]);
            title(['im = ',num2str(im),'   ','cell = ',num2str(c),'      ',data(im).file],'Interpreter','none');
            drawnow
        end
        % Lengths(i) = Length_sol;
        data(im).length_fixed(c) = Length_sol;
        data(im).width_fixed(c) = constants.width;
        if isempty( Curv_sol )
            stopodf
        end
        data(im).curvature_fixed(c) = Curv_sol;
        % pause(1)
        if donut
            %             pause
        end
        
    end
    
end
