runs = ["C:\Users\rudi\Desktop\RD\temp\E 10 eta 500 d 0.05 g 50 link_breakage dist 0.2.mat",...
"Z:\temp\E 10 eta 500 d 0.05 g 0 link_breakage dist 0.mat",...    
  "Y:\temp\E 10 eta 500 d 0.05 g 50 link_breakage dist 0.mat",...
"X:\temp\E 10 eta 500 d 0.05 g 500 link_breakage dist 0.mat",...
    ];

runs = ["C:\Users\rudi\Desktop\RD\temp\E 10 eta 500 d 0.05 g 50 link_breakage dist 0 reltol 1E-5 abstol 1E-8.mat",...
"Z:\temp\E 10 eta 500 d 0.05 g 0 link_breakage dist 0.mat",...    
  "Y:\temp\E 10 eta 500 d 0.05 g 50 link_breakage dist 0.mat",...
"X:\temp\E 10 eta 500 d 0.05 g 500 link_breakage dist 0.mat",...
];

runs = ["Z:\temp\E 10 eta 500 d 0.05 g 0 link_breakage dist 0.mat",...   
    "Y:\temp\E 10 eta 500 d 0.05 g 50 link_breakage dist 0.mat",...
  "C:\Users\rudi\Desktop\RD\temp\E 10 eta 500 d 0.05 g 50 link_breakage dist 0 reltol 1E-5 abstol 1E-8.mat",... 
  "C:\Users\rudi\Desktop\RD\temp\E 10 eta 500 d 0.05 g 50 link_breakage dist 0 reltol 5E-4 abstol 5E-4 epsilon 4E-4.mat",... 
  "Y:\temp\E 10 eta 500 d 0.05 g 50 link_breakage dist 0 reltol 5E-4 abstol 5E-4 epsilon 4E-3.mat",... 
];

names = ["g = 50 with breakage","g = 0 no breakage","g = 50 no breakage","g = 500 no breakage"];
names = ["g = 50 refined ode45","g = 0","g = 50","g = 500"];
names = ["g = 0", "g = 50", "g = 50 refined ode45", "g = 50 epsilon*10", "g = 50 epsilon*100" ];

% styles = ["-",":","-.","--"];
colors = distinguishable_colors(length(names));
tfinal = 2;
torque = 1000;
figure(438); clf;  tiledlayout(6,1);
run_inds = [1 2 4 5 ];
% run_inds = 1:length(names);
smoothing = false;
clear windows
c = 0;
for i = run_inds
   i
   
    load(runs(i),'stored_output');
    
    for j = 1:6
        
        switch j
            case 1 % speed
%                 y = sqrt(sum(stored_output.derivatives(:,1:3).^2,2));
                  y = sqrt(sum(stored_output.derivatives(:,1).^2,2));
                lims = [3.5 22];
                label = '|U| (\mum s^{-1})';
            case 2 % tail freq
                y = stored_output.derivatives(:,7);
                lims = [410 460];
                label = ('tail rad/s');
            case 3 % microns traveled / revolution
%                 y = sqrt(sum(stored_output.derivatives(:,1:3).^2,2)) ./ ( stored_output.derivatives(:,7) * 1/(2*pi) );
                 y = sqrt(sum(stored_output.derivatives(:,1).^2,2)) ./ ( stored_output.derivatives(:,7) * 1/(2*pi) );
                lims = [0.05 0.3];
                  label = ('\mum / tail rev');
%             case 4 % efficiency = speed^2 / (freq * torque)
% %                 y =  (sum(stored_output.derivatives(:,1:3).^2 , 2)) ./ (stored_output.derivatives(:,7) .* torque );
%   y =  (sum(stored_output.derivatives(:,1).^2 , 2)) ./ (stored_output.derivatives(:,7) .* torque );
%                 lims = [0 1E-3];
%                  label = '|U|^2 / (freq * torque)';
                
                  case 4 % how many network nodes are inside mesh
%                 y =  (sum(stored_output.derivatives(:,1:3).^2 , 2)) ./ (stored_output.derivatives(:,7) .* torque );
fun = @(is_inside,submesh_index) sum(submesh_index == 1 & is_inside);  % count only inside body
% fun = @(is_inside,submesh_index) sum(submesh_index == 2 & is_inside);  % count only inside tail
% fun = @(is_inside,submesh_index) sum( (submesh_index == 1 | submesh_index == 2) & is_inside);  % count inside body or tail
y = zeros(length(stored_output.time),1);  % default is none are inside
 y(stored_output.repulsion.any_flagged) = cellfun( fun, stored_output.repulsion.is_inside, stored_output.repulsion.submesh_index);
  
                 lims = [0 Inf];
                 label = '# nodes inside';
                 
            case 5 % PE
                y = stored_output.PE;
                lims = [0 2E3];
                label = 'PE';
            case 6 % repulsion force on swimmer
                y = stored_output.repulsion.total_force_mag;
                lims = [0 270];
                label = 'repulsion';
        end
        if smoothing
            [y, window] = smoothdata(y,'movmean',20);
            c = c+1;
            windows(c) = window;
        end
        
    nexttile(j)
    
% ax2 = axes('position',  [ 0.068765      0.37856     0.8875       0.1551],'FontSize',11);
% yyaxis left
% velocity_plot = plot(stored_output.time(1:1),sqrt(sum(stored_output.derivatives(1:1,1:3).^2,2)),'-','LineWidth',2);
 plot(stored_output.time,y,'LineWidth',2,'Color',colors(i,:));
% ylabel('|U| (\mum s^{-1})');
% xlabel('time (s)');
xlim([0 1.45]);
ylim(lims);
ylabel(label);
hold on
grid on

set(gca,'FontSize',12);
if j == 6 && i == run_inds(end)
    legend(names(run_inds),'Location','best','FontSize',14);
end

% ylabel('|U|^2 / power','color','r');
% ylabel('PE');
% xlabel('time (s)');  
% ylabel('repulsive force');

% repulsion_force_plot = plot(stored_output.time,stored_output.repulsion.total_torque_mag,'-','linewidth',2);
% ylabel('repulsive torque');

    end
%     pause
end
