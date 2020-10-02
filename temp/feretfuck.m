


% constants.width = 0.435;
% constants.curvature = 0.763;    % set curvature < 1E-10 to 1E-10 to avoid numerical precision problem whereby the code breaks for very straight rods
% constants.mean_width = constants.width;
min_Length = constants.mean_width;
max_Length = min( 2*pi*(1/constants.curvature + constants.mean_width/2)  ,  SF1_max*constants.mean_width); 
Lengths = linspace(min_Length, max_Length,150);
%   Lengths = 563;
Ferets = NaN(size(Lengths));  PoleDists = NaN(size(Lengths));

for L = 1:length(Lengths)
    
    if do_width_correction
        constants.width = 2*Lengths(L)/(pi - 4)*( sqrt( 1 + (pi - 4)/Lengths(L)*constants.mean_width ) - 1 );  % from mean width to actual width
    end

    
    constants.length = Lengths(L) - constants.width;
    
    if constants.length < 0
        continue
    end
    
    u0 = 0;
    [~,~,~,~,u1,u2,u3,shift, midpt] = curved_rod_parameters(constants);
    
%        u = unique( [ linspace(0,1,500)  u0 u1 u2 u3 1  ] );
         u =  linspace(0,1 - 1E-3,500)  ;   % turns out that intersections() (and InterX() ?) occasionally fails at the starting/ending pt, so just leave a tiny gap
 
    [x,y] = curved_rod_pts(u,constants,shift);
%     shat = InterX([x; y]);
    
    
%      if ~isempty(  InterX([x; y])  ) % self intersections detected
         if ~isempty(  intersections(x,y) ) % self intersections detected
%          pause
         continue
     end
    %%
    
    u_poles = [u1/2  (u2+u3)/2];  % the poles
    u_opp = u2;  % pt "opposite" u = 0 in the sense of the parameter u
    % u_half = [(u1+u2)/2  (u3+1)/2]; % lower midpt, upper midpt
    % u_guess = u_poles;
    % u_guess = [0.4 .5];
    lb = [0         u_opp];
    ub = [u_opp     1    ];
    n_guesses = 3;
    u_guesses = [linspace(lb(1),ub(1),n_guesses)' linspace(lb(2),ub(2),n_guesses)'];
    
%     options = optimoptions('fmincon', 'useparallel',false,'display','off');
    options = optimoptions('patternsearch', 'useparallel',false,'display','off');
    clear u_sol feret
    for g = 1:size(u_guesses,1)
        [u_sol(g,:),feret(g)] = patternsearch(@(X) obj_fun_feret(X,constants,shift), u_guesses(g,:), [],[],[],[],lb,ub,[],options);
    end
   
    [feret_best,ind] = min(feret);  feret_best = -feret_best;
   
    u_sol_best = u_sol(ind,:);
    [x,y] = curved_rod_pts(u_sol_best,constants,shift);
    
    %%
    u = unique( [ linspace(0,1,500)   u0 u1 u2 u3 1 ] );
    
    [x,y] = curved_rod_pts(u,constants,shift);
    
     
     
    figure(45);  plot(x,y,'-','linewidth',3);  hold on
    plot( [ x(u == u0) x(u == u1) ] ,  [ y(u == u0) y(u == u1) ] ,'k:','linewidth',1.5);
    plot( [ x(u == u2) x(u == u3) ] ,  [ y(u == u2) y(u == u3) ] ,'k:','linewidth',1.5);
    [x,y] = curved_rod_pts(u_sol_best,constants,shift);
    plot(x,y,'k--o','linewidth',1.5,'markerfacecolor','k');
    
    cx = 0 + constants.width/2;
        cy = 0 + 1/constants.curvature;
        theta_range = (Lengths(L) - constants.width) / (1/constants.curvature);
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
        
        plot(midpt(1),midpt(2),'ro','markersize',10,'markerfacecolor','r');
        
        
        
    hold off
    axis equal
    grid on
    title(['ind = ',num2str(L),'     ','feret = ',num2str(feret_best)]);
    drawnow
    
    % (1/constants.curvature + constants.width/2)*2 - feret_best
    Ferets(L) = feret_best;
    
    constants_temp = constants;  constants_temp.width = constants.mean_width;
   [PoleDists(L)] = Length2PoleDist(Lengths(L), constants_temp, false, do_width_correction);
%      pause
%     figure(35)
% plot(Lengths,PoleDists,'o-');
% grid on
% drawnow
   
end
%%
figure(34)
plot(Lengths,Ferets,'o-');
grid on

figure(35)
plot(Lengths,PoleDists,'o-');
grid on
