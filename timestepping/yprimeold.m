function [der] = yprime(t,y,Mesh,matrix_props,assembly_input,A0,A_motor_torque0,timerval,interval,skip_rigid_integrals)
%y(1:3) is current refpoint
global n_out solf
%y(4 5 6) = x, y, z rotation angles of body

%y(7) = rotation angle of tail

% body orientation described by vectors a and b, originally along the x and
% y axes at t = 0.


%% translate and rotate Mesh to current position and orientation
doplot = false;

    if length(y) == 6 %constant rotation rate condition, calculate current tail phase angle here
        y(7) =  t * assembly_input.tail.motor_freq;  %assumes initial angle = 0
    end

    y0 = [Mesh(1).refpoints(:,1); [0 0 0 0]' ];  %use initial body reference point by convention, and assume all initial angles were zero



if ~isequal(y,y0) %location and/or orientation has changed; this is not the first yprime calculation
    
      %first rotate only tail to handle tail phase angle relative to body
      Mesh(2) = rotateMesh(Mesh(2),[ y(7) 0 0]');  %rotate tail by y(7) around x-axis
%then move all submeshes according to translation and rotation of body +
%tail
    Mesh = move_Mesh(Mesh,[y(1:3) - Mesh(1).refpoints(:,1); y(4:6)]); 
    %yprime keeps track of body refpoint by convention, so define translation shift using original body refpoint
  
    
    if doplot
        fact = 1.5;
        figure(1245)
        %%
        clf
        colors = repmat([0.2 0.8 0.5],size(Mesh.elems,1),1);
        %colors(j,:) = [1 0.5 0];
        % bounds = 1E-6*[-500 500; -500 500; -500 500]*5;
        p = patch('faces',Mesh.elems(:,1:3),'vertices',Mesh.verts,'facevertexcdata',colors,'facecolor','flat','edgealpha',0.6);
        %set(p,'faceAlpha',0.6);
        axis equal
        box on
        grid on
        shading faceted
        hold on
        xlabel('x');  ylabel('y');  zlabel('z');
        % xlim([bounds(1,:)]);
        va = quiver3(y(1),y(2),y(3),a(1)*fact,a(2)*fact,a(3)*fact,'r','linewidth',3);
        vb = quiver3(y(1),y(2),y(3),b(1)*fact,b(2)*fact,b(3)*fact,'b','linewidth',3);
        vc = quiver3(tailref(1), tailref(2), tailref(3),c(1)*fact,c(2)*fact,c(3)*fact,'g','linewidth',3);
        vd = quiver3(y(1),y(2),y(3),d(1)*fact,d(2)*fact,d(3)*fact,'k','linewidth',3);
        hold off
        set(gca,'view',[-20 18]);
        drawnow
    end
    %%

    %% given current Mesh, solve for U, Omega, and possibly omega

[A,A_motor_torque] = matrix_assembly_mexed(Mesh,matrix_props,assembly_input,skip_rigid_integrals); 
            
    if skip_rigid_integrals  %rigid integrals (i.e., collocation points and elements on same rigid object) were skipped, so need to fill them in based on rotation of A0
        %A0 has to be orig t = 0 orientation that current orientations in y are in reference to
        A = rotate_A(A,y(4:6),Mesh(1),A0);  %rotate body entries
        A = rotate_A(A,[y(4)+y(7) y(5) y(6)],Mesh(2),A0); %rotate tail entries; recall that tail rotation in x is combo of body rotation and tail phase angle
    end
    
    
else %nothing at all changed, this is the first timestep
    
    A = A0;  %just use originally calculated A0, it won't change if y hasn't
    A_motor_torque = A_motor_torque0;
end

% commented this out since the current values *should* be stored in
% Mesh, after the move_Mesh operation above....
%refpoint = y(1:3);
%motor_orientation = Mesh.orientation(:,1);

[f, kinematics] = matrix_solve(Mesh, input, matrix_props, A, BCs);


solve_and_integrate_force;  %this creates f, U, Omega, and possibly omega



%%



% der(1) = dx/dt
% der(2) = dy/dt
% der(3) = dz/dt

% dx/dt = U1
% dy/dt = U2
% dz/dt = U3

der(1:3,1) = U;

der(4:6,1) = Omega;

if strcmp(motorBC,'torque')
    der(7,1) = omega;
end


if strcmp(problemtype,'freeswim')
    if strcmp(motorBC,'rotationrate')
     
        tailtorque = coeffs * f(1:n_col*3); %coeffs is output from matrix_assembly_..._mexed
        
    else
        tailtorque = motor_torque;
    end
    
end




%% save outputs besides position:  traction, speeds, and possibly tail torque



% if t >= dt_out * (n_out - 1)

% assignin('base','ftemp',f);
%assignin('base','n_out',n_out
%  evalin('base','solf(n_out,:) = ftemp(1:n_col*3);');
%solf(n_out,:) = [t; f(1:n_col*3)];

if n_out > size(solf,1)
    solf(end+1:end+50,:) = NaN;
end

switch motorBC
    case 'torque'
        solf(n_out,:) = [t; f; tailtorque]; %has rotationrate in f
    case 'rotationrate'
        solf(n_out,:) = [t; f; motor_freq; tailtorque];
end
elapsed = toc(timerval);
fraction_done = (t - interval(1)) /  (interval(2) - interval(1));
fraction_left = 1 - fraction_done;
timeleft = elapsed / fraction_done * fraction_left;
disp(['n_out = ',num2str(n_out),'       t = ',num2str(t),'       elapsed = ',num2str(elapsed/60,2),' min','       timeleft = ',num2str(timeleft/60,2),' min']);
n_out = n_out + 1;

drawnow

% end




% da/dt = W cross a

% der(4) = Omega(2)*a(3) - Omega(3)*a(2);
% der(5) = Omega(3)*a(1) - Omega(1)*a(3);
% der(6) = Omega(1)*a(2) - Omega(2)*a(1);
%
% % db/dt = W cross b
%
% der(7) = Omega(2)*b(3) - Omega(3)*b(2);
% der(8) = Omega(3)*b(1) - Omega(1)*b(3);
% der(9) = Omega(1)*b(2) - Omega(2)*b(1);
%
% %if strcmp(motorBC,'torque')
% %     omega = varargin{3};
%
%
% % dc/dt = omega cross c
%
% der(10) = -omega*c(2);
% der(11) = omega*c(1);

%end