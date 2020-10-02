namelist = {'ubercoarse','coarsest','coarser','base','finer','finest'};
%folder = 'C:\Users\rudi\Desktop\RD\code\dumps\';
folder = 'E:\Hull\code\dumps\';

parameters_refinement = 'coarser';

body_refinements = {'coarse','coarse', 'base','base',     'fine','fine'};
tail_refinements = {'base','fine_ends','base','fine_ends','base','fine_ends'  };

clear fcoeffs speeds

%for f = 1:length(namelist)
    for f = 1:length(body_refinements)
    f
   %temp = dir([folder,'*forced*_',namelist{f},'_dump.mat']);
    
    temp = dir([folder,'*0.85_',body_refinements{f},'*1.49_',tail_refinements{f},'*forced*',parameters_refinement,'_dump.mat']);
    
    
    temp = temp.name;
    
    temp = load(temp);
    try
  %temp.solutions.(flowcase).forces.drag
    %temp.input.accuracy.integration_tol.traction
    end
    fcoeffs(f) = temp.fcoeffs;
    
    
   %   temp = dir([folder,'*_',namelist{f},'_timestepping.mat']);
      
       temp = dir([folder,'*0.85_',body_refinements{f},'*1.49_',tail_refinements{f},'*motor*',parameters_refinement,'_timestepping.mat']);
       
       
    temp = temp.name;
    
    temp = load(temp);
    
      speeds(f) = temp.timestepping_convergence.avg_speed;
      
end
%%


%clear fcoeffs speeds
folder = input.paths.dumpfolder;

% namelist = { 'curved_rod_AR1_8_AR2_0.85_fine_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_base_forced_eps_8e-05_abstol_0.004_interptol_0.08_ubercoarse-refined-integration_dump',...
%     'curved_rod_AR1_8_AR2_0.85_base_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_base_forced_eps_4e-05_abstol_0.008_interptol_0.08_ubercoarse-refined-eps_dump',...
% 'curved_rod_AR1_8_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_8e-05_abstol_0.008_interptol_0.08_ubercoarse_dump',...
% 'curved_rod_AR1_8_AR2_0.85_base_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_base_forced_eps_4e-05_abstol_0.004_interptol_0.04_coarser2_dump',...
% 'curved_rod_AR1_8_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_2.5e-06_abstol_0.00025_interptol_0.0025_finest_dump' };
% 
% 
% namelist2 = {'curved_rod_AR1_8_AR2_0.85_fine_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_base_motorBC_torque_eps_8e-05_abstol_0.004_interptol_0.08_ubercoarse-refined-integration_timestepping',...
%     'curved_rod_AR1_8_AR2_0.85_base_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_base_motorBC_torque_eps_4e-05_abstol_0.008_interptol_0.08_ubercoarse-refined-eps_timestepping',...
%     'curved_rod_AR1_8_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_8e-05_abstol_0.008_interptol_0.08_ubercoarse_timestepping',...
%     'curved_rod_AR1_8_AR2_0.85_base_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_base_motorBC_torque_eps_4e-05_abstol_0.004_interptol_0.04_coarser2_timestepping',...
%     'curved_rod_AR1_8_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_2.5e-06_abstol_0.00025_interptol_0.0025_finest_timestepping' };

namelist = {...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01force_torque_0.1_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01timestepping_10_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.5interp_50_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_2e-05_abstol_0.002_interptol_0.01epsilon_2_traction_2_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_5e-06_abstol_0.0005_interptol_0.01epsilon_0.5_traction_0.5_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01timestepping_reltol_100_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.2interp_20_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_2e-05_abstol_0.001_interptol_0.01epsilon_2_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.002_interptol_0.01traction_2_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.001_interptol_0.01epsilon_4_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.004_interptol_0.01traction_4_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01timestepping_reltol_10000_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.008_interptol_0.01traction_8_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.016_interptol_0.01traction_16_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.032_interptol_0.01traction_32_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_10interp_1000_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_8e-05_abstol_0.001_interptol_0.01epsilon_8_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.1epsilon_4_traction_16_interp_10_ode45_100_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01timestepping_reltol_1000_dump',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.1epsilon_4_traction_16_interp_10_ode45_1000_dump',...  % 21
   ...
   'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01_dump',...  %22
    'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.1epsilon_4_traction_16_interp_10_ode45_100_dump',... %23
    'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.01epsilon_4_traction_16_interp_10_ode45_100_interrogation_2_dump',...  %24
    'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.01epsilon_4_traction_16_interp_10_ode45_10_interrogation_4.5_dump',... %25
    'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.01epsilon_4_traction_16_interp_10_ode45_100_interrogation_4.5_dump',... %26
    'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.01epsilon_4_traction_16_interp_10_ode45_1000_interrogation_4.5_dump',... %27
   'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.01epsilon_4_traction_16_interp_10_ode45_1_interrogation_4.5_dump',...
     'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01interrogation_4.5_dump',... % 29
   ...
   'curved_rod_AR1_12_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01_ode45_10_dump',...
   'curved_rod_AR1_12_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.1_epsilon_4_traction_16_interp_10_interrogation_4_ode45_10_dump',...
  ...
   'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01nontraction_0.1_interrogation_4.5_dump',...  %32
   'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01nontraction_0.1_interrogation_4.5_dump',... %33
   'curved_rod_AR1_12_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_1e-05_abstol_0.001_interptol_0.01nontraction_0.1_interrogation_4.5_dump',... %34
   'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.1eps_4_traction_16_interp_10_timestep_2_ode45_10_interrogation_2.25_dump',... %35
   'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.1eps_4_traction_16_interp_10_timestep_2_ode45_10_interrogation_2.25_dump',... %36
   'curved_rod_AR1_12_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_forced_eps_4e-05_abstol_0.016_interptol_0.1eps_4_traction_16_interp_10_timestep_2_ode45_10_interrogation_2.25_dump',... %37
    };
    
    
    namelist2 = {...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01force_torque_0.1_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01timestepping_10_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.5interp_50_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_2e-05_abstol_0.002_interptol_0.01epsilon_2_traction_2_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_5e-06_abstol_0.0005_interptol_0.01epsilon_0.5_traction_0.5_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01timestepping_reltol_100_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.2interp_20_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_2e-05_abstol_0.001_interptol_0.01epsilon_2_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.002_interptol_0.01traction_2_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.001_interptol_0.01epsilon_4_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.004_interptol_0.01traction_4_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01timestepping_reltol_10000_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.008_interptol_0.01traction_8_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.016_interptol_0.01traction_16_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.032_interptol_0.01traction_32_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_10interp_1000_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_8e-05_abstol_0.001_interptol_0.01epsilon_8_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.1epsilon_4_traction_16_interp_10_ode45_100_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01timestepping_reltol_1000_timestepping',...
    'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.1epsilon_4_traction_16_interp_10_ode45_1000_timestepping',...
 ...
  'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01_timestepping',...  
 'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.1epsilon_4_traction_16_interp_10_ode45_100_timestepping',...
    'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.01epsilon_4_traction_16_interp_10_ode45_100_interrogation_2_timestepping',...
    'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.01epsilon_4_traction_16_interp_10_ode45_10_interrogation_4.5_timestepping',...
    'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.01epsilon_4_traction_16_interp_10_ode45_100_interrogation_4.5_timestepping',...
    'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.01epsilon_4_traction_16_interp_10_ode45_1000_interrogation_4.5_timestepping',...
  'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.01epsilon_4_traction_16_interp_10_ode45_1_interrogation_4.5_timestepping',...
  'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01interrogation_4.5_timestepping',...
  ...
   'curved_rod_AR1_12_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01_ode45_10_timestepping',...
   'curved_rod_AR1_12_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.1_epsilon_4_traction_16_interp_10_interrogation_4_ode45_10_timestepping',...
  ...
   'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01nontraction_0.1_interrogation_4.5_timestepping',... 
   'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01nontraction_0.1_interrogation_4.5_timestepping',...
   'curved_rod_AR1_12_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_1e-05_abstol_0.001_interptol_0.01nontraction_0.1_interrogation_4.5_timestepping',...
   'curved_rod_AR1_2_AR2_0.05_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.1eps_4_traction_16_interp_10_timestep_2_ode45_10_interrogation_2.25_timestepping',...
   'curved_rod_AR1_4_AR2_0.5_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.1eps_4_traction_16_interp_10_timestep_2_ode45_10_interrogation_2.25_timestepping',...
   'curved_rod_AR1_12_AR2_0.85_tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49_motorBC_torque_eps_4e-05_abstol_0.016_interptol_0.1eps_4_traction_16_interp_10_timestep_2_ode45_10_interrogation_2.25_timestepping',...
   };
    

    
    
for f = 30:length(namelist)
  %  for f = 1:length(body_refinements)
    f/length(namelist)
   %temp = dir([folder,'*forced*_',namelist{f},'_dump.mat']);
    
  %  temp = dir([folder,'*0.85_',body_refinements{f},'*1.49_',tail_refinements{f},'*forced*',parameters_refinement,'_dump.mat']);

  name = [folder,namelist{f},'.mat'];
  
    temp = load(name);
    
    try
  %temp.solutions.(flowcase).forces.drag
    %temp.input.accuracy.integration_tol.traction
    end
    fcoeffs(f) = temp.fcoeffs;
    
    
%      temp = dir([folder,'*_',namelist{f},'_timestepping.mat']);
%       
%       % temp = dir([folder,'*0.85_',body_refinements{f},'*1.49_',tail_refinements{f},'*motor*',parameters_refinement,'_timestepping.mat']);
%        
%        
%     temp = temp.name;
    
      name = [folder,namelist2{f},'.mat'];
  
    temp = load(name);
    
      speeds(f) = temp.timestepping_convergence.avg_speed;
      
    end

    %%
    
    %for when multiple lines are reported, these are the different geom
    %cases.  When only one line reported, it is for the first case.
    % AR1 = 4   AR2 = 0.5
    %       2         0.05
    %       5         0.95

     'effect of refining force and torque abstol by 1/10'
    (fcoeffs(1).rotation - fcoeffs(2).rotation) ./fcoeffs(2).rotation *100
    (fcoeffs(1).translation - fcoeffs(2).translation) ./fcoeffs(2).translation *100
    (speeds(1) - speeds(2)) / speeds(2) * 100
    % almost no effect
    
    'effect of unrefining timestepping convergence tol by * 10'
        (fcoeffs(1).rotation - fcoeffs(3).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(3).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(3)) / speeds(1) * 100
    % 0.02% change                                          so leave as-is
    
           'effect of unrefining interp tol by * 20'
            (fcoeffs(1).rotation - fcoeffs(8).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(8).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(8)) / speeds(1) * 100
    % almost no effect
    
    'effect of unrefining interp tol by * 50'
            (fcoeffs(1).rotation - fcoeffs(4).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(4).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(4)) / speeds(1) * 100
    % almost no effect  this is a little under an interp convergence reltol
    % of 0.1, so I guess that's good enough (it yields 9 phase angle eval points)
    
              'effect of unrefining interp tol by * 1000'
            (fcoeffs(1).rotation - fcoeffs(17).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(17).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(17)) / speeds(1) * 100
    %  almost no effect  but this is just so ridiculously coarse (5 phase
    %  angle eval points) that we should probably do more
    
    'effect of unrefining epsilon and traction by * 2'
            (fcoeffs(1).rotation - fcoeffs(5).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(5).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(5)) / speeds(1) * 100
    % 0.002% change
        
    'effect of refining epsilon and traction by / 2'
            (fcoeffs(1).rotation - fcoeffs(6).rotation) ./fcoeffs(6).rotation *100
    (fcoeffs(1).translation - fcoeffs(6).translation) ./fcoeffs(6).translation *100
    (speeds(1) - speeds(6)) / speeds(6) * 100
    %   0.001% change
    
        'effect of unrefining ode45 reltol by * 100'
            (fcoeffs(1).rotation - fcoeffs(7).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(7).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(7)) / speeds(1) * 100
    %  almost no effect
    
     'effect of unrefining ode45 reltol by * 1000'
            (fcoeffs(1).rotation - fcoeffs(20).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(20).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(20)) / speeds(1) * 100
    %  almost no effect      so use * 1000 since * 10000
    %  is no good
    
            'effect of unrefining ode45 reltol by * 10000'
            (fcoeffs(1).rotation - fcoeffs(13).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(13).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(13)) / speeds(1) * 100
    %  1% error
    
    
        'effect of unrefining epsilon by * 2'
            (fcoeffs(1).rotation - fcoeffs(9).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(9).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(9)) / speeds(1) * 100
    % 0.003% change
    
                'effect of unrefining epsilon by * 4'
            (fcoeffs(1).rotation - fcoeffs(11).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(11).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(11)) / speeds(1) * 100
    % 0.008% change
    
                'effect of unrefining epsilon by * 8'
            (fcoeffs(1).rotation - fcoeffs(18).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(18).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(18)) / speeds(1) * 100
    % 0.02% change
    
    
            'effect of unrefining traction by * 2'
            (fcoeffs(1).rotation - fcoeffs(10).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(10).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(10)) / speeds(1) * 100
    % almost no change
    

    
            'effect of unrefining traction by * 4'
            (fcoeffs(1).rotation - fcoeffs(12).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(12).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(12)) / speeds(1) * 100
    % almost no change
    
               'effect of unrefining traction by * 8'
            (fcoeffs(1).rotation - fcoeffs(14).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(14).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(14)) / speeds(1) * 100
    % 0.003% change
    
               'effect of unrefining traction by * 16'
            (fcoeffs(1).rotation - fcoeffs(15).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(15).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(15)) / speeds(1) * 100
    % 0.01% change
    
               'effect of unrefining traction by * 32'
            (fcoeffs(1).rotation - fcoeffs(16).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(16).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(16)) / speeds(1) * 100
    % 0.03% change
    
                 'effect of unrefining epsilon by * 4, traction * 16, interp * 10, ode45 * 100'
            (fcoeffs(1).rotation - fcoeffs(19).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(19).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(19)) / speeds(1) * 100
    % 0.01% change     woot!       takes 6 minutes, timestepping takes 8 min
     % 0.01% change     woot!       takes 6 minutes, timestepping takes 23 min
     %      change                  takes 7 minutes, timestepping takes 11 min
    
                     'effect of unrefining epsilon by * 4, traction * 16, interp * 10, ode45 * 1000'
            (fcoeffs(1).rotation - fcoeffs(21).rotation) ./fcoeffs(1).rotation *100
    (fcoeffs(1).translation - fcoeffs(21).translation) ./fcoeffs(1).translation *100
    (speeds(1) - speeds(21)) / speeds(1) * 100
    % 0.01% change     woot!    takes 6 minutes, timesteppping takes 4.4 min
    
    
    
    
    %%  AR1 = 2   AR2 = 0.05
    
     'effect of unrefining epsilon by * 4, traction * 16, interp * 10, ode45 * 100'
            (fcoeffs(22).rotation - fcoeffs(23).rotation) ./fcoeffs(22).rotation *100
    (fcoeffs(22).translation - fcoeffs(23).translation) ./fcoeffs(22).translation *100
    (speeds(22) - speeds(23)) / speeds(22) * 100
    % 0.02% change     
    
    
       'effect of interrogation time * 2'
            (fcoeffs(23).rotation - fcoeffs(24).rotation) ./fcoeffs(23).rotation *100
    (fcoeffs(23).translation - fcoeffs(24).translation) ./fcoeffs(23).translation *100
    (speeds(23) - speeds(24)) / speeds(23) * 100
    % 0.002% change     
    
       'effect of refining ode45 reltol / 10 and refine interrogation time * > 2'
            (fcoeffs(25).rotation - fcoeffs(24).rotation) ./fcoeffs(25).rotation *100
    (fcoeffs(25).translation - fcoeffs(24).translation) ./fcoeffs(25).translation *100
    (speeds(25) - speeds(24)) / speeds(25) * 100
    % almost no change
    
      ' refine interrogation time * > 2'
            (fcoeffs(24).rotation - fcoeffs(26).rotation) ./fcoeffs(26).rotation *100
    (fcoeffs(24).translation - fcoeffs(26).translation) ./fcoeffs(26).translation *100
    (speeds(24) - speeds(26)) / speeds(26) * 100
    % almost no change

             ' unrefine ode45 reltol by * 10 (* 10 to * 100)'
            (fcoeffs(26).rotation - fcoeffs(25).rotation) ./fcoeffs(25).rotation *100
    (fcoeffs(26).translation - fcoeffs(25).translation) ./fcoeffs(25).translation *100
    (speeds(26) - speeds(25)) / speeds(25) * 100
    % almost no change
    
         ' unrefine ode45 reltol by * 10 (* 100 to * 1000)'
            (fcoeffs(26).rotation - fcoeffs(27).rotation) ./fcoeffs(26).rotation *100
    (fcoeffs(26).translation - fcoeffs(27).translation) ./fcoeffs(26).translation *100
    (speeds(26) - speeds(27)) / speeds(26) * 100
    % 0.1% change
    
    
    'unrefine epsilon * 4, traction * 16, interp tol * 10.   interrogation * 4.5.  ode45 * 1'
          (fcoeffs(29).rotation - fcoeffs(28).rotation) ./fcoeffs(29).rotation *100
    (fcoeffs(29).translation - fcoeffs(28).translation) ./fcoeffs(29).translation *100
    (speeds(29) - speeds(28)) / speeds(29) * 100
    
    
     'unrefine epsilon * 4, traction * 16, interp tol * 10.   interrogation * 4.5.  ode45 * 10'
          (fcoeffs(29).rotation - fcoeffs(25).rotation) ./fcoeffs(29).rotation *100
    (fcoeffs(29).translation - fcoeffs(25).translation) ./fcoeffs(29).translation *100
    (speeds(29) - speeds(25)) / speeds(29) * 100
    
    
     'unrefine epsilon * 4, traction * 16, interp tol * 10.   interrogation * 4.5.  ode45 * 100'
          (fcoeffs(29).rotation - fcoeffs(26).rotation) ./fcoeffs(29).rotation *100
    (fcoeffs(29).translation - fcoeffs(26).translation) ./fcoeffs(29).translation *100
    (speeds(29) - speeds(26)) / speeds(29) * 100
    
     'unrefine epsilon * 4, traction * 16, interp tol * 10.   interrogation * 4.5.  ode45 * 1000'
          (fcoeffs(29).rotation - fcoeffs(27).rotation) ./fcoeffs(29).rotation *100
    (fcoeffs(29).translation - fcoeffs(27).translation) ./fcoeffs(29).translation *100
    (speeds(29) - speeds(27)) / speeds(29) * 100
    
    %% AR1 = 12   AR2 = 0.85
    
    %lack of convergence of avg speed unless ode45 reltol < * 100
    %unrefinement. ode45 * 10 and * 1 look the same, ~ no change in avg_speed
    
       'epsilon * 4, traction * 16, interp * 10,  interrogation  / 4'
          (fcoeffs(31).rotation - fcoeffs(30).rotation) ./fcoeffs(30).rotation *100
    (fcoeffs(31).translation - fcoeffs(30).translation) ./fcoeffs(30).translation *100
    (speeds(31) - speeds(30)) / speeds(30) * 100
    % 0.01 % change
    %% Final comparisons:  each has 0.01 - 0.02% max change  WOOT

% AR1 = 2  AR2 = 0.05
       'epsilon * 4, traction * 16, other integration * 10, interp * 10, timestep conv * 2 ode45 *10 interrogation decreased from * 4.5 to * 2.25'
          (fcoeffs(32).rotation - fcoeffs(35).rotation) ./fcoeffs(32).rotation *100
    (fcoeffs(32).translation - fcoeffs(35).translation) ./fcoeffs(32).translation *100
    (speeds(32) - speeds(35)) / speeds(32) * 100
    % 0.01 % change
    
    % AR1 = 4  AR2 = 0.5
       'epsilon * 4, traction * 16, other integration * 10, interp * 10, timestep conv * 2 ode45 *10 interrogation decreased from * 4.5 to * 2.25'
          (fcoeffs(33).rotation - fcoeffs(36).rotation) ./fcoeffs(33).rotation *100
    (fcoeffs(33).translation - fcoeffs(36).translation) ./fcoeffs(33).translation *100
    (speeds(33) - speeds(36)) / speeds(33) * 100
    % 0.01 % change
    
        % AR1 = 12  AR2 = 0.85
       'epsilon * 4, traction * 16, other integration * 10, interp * 10, timestep conv * 2 ode45 *10 interrogation decreased from * 4.5 to * 2.25'
          (fcoeffs(34).rotation - fcoeffs(37).rotation) ./fcoeffs(34).rotation *100
    (fcoeffs(34).translation - fcoeffs(37).translation) ./fcoeffs(34).translation *100
    (speeds(34) - speeds(37)) / speeds(34) * 100
    % 0.01 % change
    %%

diffs.speed = abs((speeds - speeds(end)) ./ speeds ) * 100;
diffs.translation = abs( (vertcat(fcoeffs.translation) - repmat(fcoeffs(end).translation,length(fcoeffs),1)) ./ vertcat(fcoeffs.translation) ) * 100;
diffs.rotation = abs( (vertcat(fcoeffs.rotation) - repmat(fcoeffs(end).rotation,length(fcoeffs),1)) ./ vertcat(fcoeffs.rotation) ) * 100;

diffs.speed
diffs.translation
diffs.rotation

