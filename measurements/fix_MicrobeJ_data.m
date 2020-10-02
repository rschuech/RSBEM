% widths = test_data.width;
% curvatures = test_data.curvature;
% Ferets = test_data.length;
% Lengths = NaN(size(widths));

% fix_length_bug = false;  % true for old MicrobeJ versions with Feret -
% Length bug, false for new MicrobeJ that outputs correct Length - gets set
% automatically later
do_width_correction = true;  % convert measured mean width to geometric curved rod width
do_curv_correction = true;  % convert measured curvature2 based on pole points to correct curvature
% do_width_correction = false;


% redo_files = {'Acidimicrobium ferrooxidans gen nov sp nov mixed-culture ferrous iron oxidation with Sulfobacillus species_Page_3_Image_0003'};
% redo_files = {'Desulfobacter vibrioformis - Lien and Beeder 1997_Page_2_Image_0002_edited',...
%     'Desulfonatronum zhilinae sp nov a novel haloalkaliphilic sulfate-reducing bacterium from soda Lake Alginskoe Trans-Baikal Region Russia_Page_4_Image_0001 panel B',...
%     'Evaluation of Leptospirillum ferrooxidans for leaching_Page_4_Image_0002 panel A',...
%     'Interfacial activity and leaching patterns of Leptospirillum ferrooxidans on pyrite_Page_03_Image_0001',...
%     'Leptospirillum ferriphilum - Liu et al 2007_Page_4_Image_0002.txt',...
%     'New types of acetate-oxidizing sulfate-reducing Desulfobacter species D hydrogenophilus sp nov D latus sp nov and D curvatus sp nov_Page_3_Image_0002_Page_1_Image_0001 panel D',...
%     'vol2C-nnn-910-1114'};

% redo_files = {'vol2C-nnn-910-1114',...
%     'Shewanella halifaxensis - Zhao et al 2006_Page_2_Image_0001_B_edited',...
%     'Leptospirillum ferrooxidans_Giavenoetal2007'};

do_plots = false;
% do_plots = true;


% bad_inds = [];

ppm = ParforProgressStarter2('fixing MicrobeJ results', length(data));

for im = 1:length(data) % images
    %     data(im).file
    
    Length_sol = []; max_Length = [];  max_Feret = [];
    
    if ~isfield(data(im),'length_fixed') || isempty(data(im).length_fixed)
        data(im).length_fixed = NaN(size(data(im).length));
        data(im).width_fixed = NaN(size(data(im).width));
        data(im).curvature_fixed = NaN(size(data(im).curvature));
        data(im).donut = NaN(size(data(im).length));
        data(im).fixed_yet = false(size(data(im).length));
        data(im).SF1 = false(size(data(im).length));
        data(im).SF2 = false(size(data(im).length));
    end
    
    

    
    %     donuts = false(1,length(data(im).length));
    Curv_sol = []; u_sol_best = [];
    for c = 1:length(data(im).length)  %individual cells
        %         c
        
        
        
        
        if ~ismember(data(im).file,redo_files)  && ( data(im).fixed_yet(c) )
            continue
        end
        
        if data(im).feret_fixed(c)
            fix_length_bug = false;
        else
            fix_length_bug = true;
        end
        
        
        
        
        donut_safety_factor = 0.05;  % if known Feret value is within this factor of the Feret of a semicircular-or-greater rod, then switch to PoleDist method
        SF1_max = 50; % no way SF1 can ever be larger than this
        
        constants = [];
        constants.width = data(im).width(c);
        
        
        constants.curvature = data(im).curvature(c);
        
        Feret = data(im).length(c);
        %         PoleDist = data(im).poledist(c);
        
        %         PoleDist = NaN;
        
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
        
        
        
        
%         if fix_length_bug
            objfun = @(Length) -Length;
            %         Length_guess = (global_min_Length + global_max_Length) / 2;
            Length_guess = global_min_Length * 1;
            % max_Length = fmincon(objfun,Length_guess, [],[],[],[],global_min_Length, global_max_Length, nonlincon,options);
            [max_Length,~,flag] = patternsearch(objfun,Length_guess, [],[],[],[],global_min_Length, global_max_Length, nonlincon,options);
            if flag <= 0
                error('max_Length calc failed');
            end
            max_Feret = Length2Feret(max_Length, constants,true, do_width_correction);
%         end
        %%
        
        if Feret < (1 - donut_safety_factor) * max_Feret
            data(im).donut(c) = false;
        else
            data(im).donut(c) = true;
        end
        
        
        if fix_length_bug
            if  ~data(im).donut(c)
                Feret_Curv2 = [Feret , constants.curvature];
                check_intersections = true;  % makes almost no difference in speed so safer to do it
                objfun = @(Length_Curv) objfun_feret_curv(Length_Curv, Feret_Curv2, constants, check_intersections, do_width_correction);
                lb = [min_Length , constants.curvature * 1E-1];  ub = [max_Length , constants.curvature * 3];
                %             guess = [Length_sol, constants.curvature];
                guess = [Feret, constants.curvature* 1];  % * 3
                %                 guess = [258, constants.curvature* 1.2];  % * 3
                
                %                 options = optimoptions('patternsearch','StepTolerance',1E-6/1,'MeshTolerance',1E-6/1,'FunctionTolerance',1E-12,'MaxIterations',1E6,'ConstraintTolerance',1E-8,'display','iter');
                
                %                 [Length_Curv,fval,flag] = patternsearch(objfun,guess, [],[],[],[],lb, ub, [],options);
                
                options = optimoptions('fsolve','FunctionTolerance',1E-10,'OptimalityTolerance',1E-10, 'MaxFunctionEvaluations', 1E5, 'MaxIterations',1E4,'display','off');
                [Length_Curv,fval,flag] = fsolve(objfun, guess,options);
                
                if any((abs( Length_Curv - lb ) ./ lb ) < 0.01) || any((abs(ub - Length_Curv) ./ ub) < 0.01)
                    disp(['im = ',num2str(im),'   ','c = ',num2str(c)]);
                    
                    error('Length and Curvature solution too close to bounds');
                end
                
                if flag <= 0
                    error('Length and Curvature fix calc completely failed');
                end
                
                [~, Feret_sol, Curv2_sol] = objfun(Length_Curv);
                ansy = [Feret_Curv2;  Feret_sol Curv2_sol];
                
                abs_errors = abs( [Feret - Feret_sol , constants.curvature - Curv2_sol ] );
                
                if any(abs_errors > 1E-9)  % was 1E-6
                    error('Abs errors too big');
                end
                
                Length_sol = Length_Curv(1);
                Curv_sol = Length_Curv(2);
                
            else
                
                %             donuts(c) = true;
                %             disp(['Feret cutoff exceeded for ',data(im).file,'    ','cell index ',num2str(c)]);
                %             bad_inds = unique([bad_inds   im]);
                
                %             if isnan(PoleDist)
                Length_sol = NaN;
                Curv_sol = NaN;
                disp(['Donut:  ',data(im).file,'    ','cell index ',num2str(c)]);
                data(im).length_fixed(c) = Length_sol;
                data(im).width_fixed(c) = NaN;
                data(im).curvature_fixed(c) = Curv_sol;
                data(im).fixed_yet(c) = false;
                
                continue
            end
            
            
        else  % no need to fix length bug but still fix width and curvature
            %                     [obj, Curvature2] = objfun_curv(Curv, Curv2, constants, check_intersections, do_width_correction);
            constants.Length = data(im).length(c);
            
            Curv2 =  constants.curvature;
            check_intersections = true;  % makes almost no difference in speed so safer to do it
            objfun = @(Curv) objfun_curv(Curv, Curv2, constants, check_intersections, do_width_correction);
            lb = [ constants.curvature * 1E-1];  ub = [ constants.curvature * 3];
            %             guess = [Length_sol, constants.curvature];
            guess = [ constants.curvature* 1];  % * 3
            %                 guess = [258, constants.curvature* 1.2];  % * 3
            
            %                 options = optimoptions('patternsearch','StepTolerance',1E-6/1,'MeshTolerance',1E-6/1,'FunctionTolerance',1E-12,'MaxIterations',1E6,'ConstraintTolerance',1E-8,'display','iter');
            
            %                 [Length_Curv,fval,flag] = patternsearch(objfun,guess, [],[],[],[],lb, ub, [],options);
            
            options = optimoptions('fsolve','FunctionTolerance',1E-10,'OptimalityTolerance',1E-10, 'MaxFunctionEvaluations', 1E5, 'MaxIterations',1E4,'display','iter');
            [Curv_sol,fval,flag] = fsolve(objfun, guess,options);
            Length_sol = constants.Length;
            
            if flag <= 0
                error('Curvature fix calc completely failed');
            end
            
            if any((abs( Curv_sol - lb ) ./ lb ) < 0.05) || any((abs(ub - Curv_sol) ./ ub) < 0.05)
                error('Length and Curvature solution too close to bounds');
            end
            
            
            [~, Curv2_sol] = objfun(Curv_sol);
            abs_error = abs(Curv2_sol - Curv2);
            if (abs_error > 1E-9)
                error('Abs error too big');
            end
            
            
            
            
            
            
            width = 2*Length_sol/(pi - 4)*( sqrt( 1 + (pi - 4)/Length_sol*constants.width ) - 1 );
            
            shaperoo = [  Length_sol / width   ,     Length_sol / (2*pi*1/Curv_sol) ]
            
            % computed minlength 131.93    fval 7.86E-18    length
            % 133.87
            
            % computed minlength 103.62    fval 3.55E-18    length
            % 127.48
            
        end
        constants2 = constants;  constants2.curvature = Curv_sol;
        [~,u_sol_best] = Length2Feret(Length_sol, constants2, true, do_width_correction);
        
        
        
        
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
        
        
        
        data(im).length_fixed(c) = Length_sol;
        data(im).width_fixed(c) = constants.width;
        if isempty( Curv_sol )
            stopodf
        end
        data(im).curvature_fixed(c) = Curv_sol;
        data(im).fixed_yet(c) = true;
        % pause(1)
        %         if donut
        %             pause
        %         end
        
    end
    
        data(im).SF1 = data(im).length_fixed ./ data(im).width_fixed;
    data(im).SF2 = data(im).length_fixed ./ (2*pi* 1./data(im).curvature_fixed);


    ppm.increment(im);
end

delete(ppm);
