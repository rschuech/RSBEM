% shat 2 get working as mex:  save_Col_inds_mexed , store mesh constants,
% matrix assembly mexed, compute_BCs_mexed

% also, investigate 2% difference in field_vel from
% field_velocity.m vs field_velocity_mexed

% improve T6interp by using simplified formulas if all alpha, beta, gamma
% for the Mesh = 0.5 (precompute this once beforehand and send in a flag).

% check if shat is converged vs integration tol for field velocity for dino

%ALL SETTINGS AND INPUT PARAMETERS

%any values that are in column vector N x 1 {cell arrays} are incorporated into a parameter
%sweep that covers all combinations of each swept parameter.

%any values that are in row vector 1 X N {cell arrays} are not swept
%through all possible combinations, but simply iterated through one by one
%this means that all 1 x N cell array inputs MUST be the same length, since
%they must all match up

%Values not in {cell arrays}, whether they are scalar, vector, or matrix, are assumed to be
%individual (possibly non-scalar) parameters (e.g. a 3 x 1 orientation vector) and are kept as
%single parameter values throughout

%inputs is what is set directly below, with cell arrays for swept values
%Inputs is created at the end, with one combination of all swept values per
%index

clear inputs Inputs

switch getenv('computername')
    case 'UBERTOP'
        % cd E:/Hull/git_code/;
        cd C:\Hull\git_code2\;
        inputs.paths.datfolder = 'C:\Hull\body meshes\';   %location of .dat T6 mesh file(s)
       % inputs.paths.intersections_file = 'E:/Hull/git_code/submesh_intersections.mat';
        rack = getenv('computername');
        inputs.paths.rack = rack;
         inputs.paths.results_folder = 'C:\Hull\temp\';
         inputs.paths.results_file = [rack,'.mat'];
        inputs.paths.dumpfolder = 'C:\Hull\temp\';
        
    case {'CFD01','CFD02','CFD03','CFD04'}
        cd C:\Users\rudi\Desktop\RD\git_code2\;
        
        
        inputs.paths.datfolder = 'C:\Users\rudi\Desktop\RD\dino_mesh_automation\phases_hairs_2_1.5\final\';
        prefix = '';
        
%         inputs.paths.intersections_file = 'C:\Users\rudi\Desktop\RD\submesh_intersections.mat';
        rack = getenv('computername');
        inputs.paths.rack = rack;
        inputs.paths.results_folder = 'C:\Users\rudi\Desktop\RD\pape main results\coplanar_2_again\';
        %         inputs.paths.results_file = results_file_name;
        
        inputs.paths.results_file = [rack,'.mat'];
        
        inputs.sync_results = false;  %check all Results files for current geometry and overwrite / append if found (if false, only do this for local Results file)
        inputs.paths.lock_file = 'dino_lock'; % to avoid race condition between racks when writing to shared Results file
        if strcmp(getenv('computername'), 'CFD01')
            inputs.paths.global_lock_file = 'C:\Users\rudi\Desktop\RD\Results\global_lock';
        else
            inputs.paths.global_lock_file = 'X:\Results\global_lock';
        end
        distcomp.feature( 'LocalUseMpiexec', false ); %fixes bug in parallel toolbox as of 2015a
        %         inputs.paths.dumpfolder = 'C:\Users\rudi\Desktop\RD\cleaner_dumps\';
        %          inputs.paths.dumpfolder = 'C:\Users\rudi\Desktop\RD\dino_dumps\';
        inputs.paths.dumpfolder = inputs.paths.results_folder;
%          inputs.paths.dumpfolder = 'C:\Users\rudi\Desktop\RD\dino_dumps_hairs_2_1.5_no_offset\';
        
end




digits = 16;

inputs.use_interp_partial = true;
inputs.interp_phases_tol = -2;  %round phase values to 10^(this) when looking for matches in interp_partial file etc
% inputs.discard_interp_phases = [3.338]; %whether or not we did these before, discard them and don't redo; it isn't allowed to discard 0 since this is where the periodicity is enforced
% inputs.discard_interp_phases = [5.792 2.651 0.2945];
inputs.discard_interp_phases = [];
% inputs.redo_interp_phases = [0 0.098 0.196 0.2945]; % if any of these are found in interp_partial file, redo them instead of using old data
inputs.redo_interp_phases = [];
inputs.plot_discarded_interp_phases = false; % plot old stored values that are now discarded in spline fits?
inputs.kinematics_interp_method = 'spline';  % either 'trig' for Fourier interpolation (won't work if any interp phases are discarded!) or 'spline' (still works with discarded points)
inputs.periodic_kinematics = true;
assert(all(~ismember(roundn([0 2*pi],inputs.interp_phases_tol),roundn(inputs.discard_interp_phases,inputs.interp_phases_tol))) , 'Not allowed to discard phase 0 or 2*pi since this is where periodicity is enforced');

addpath(genpath('./')); %add all subfolders to path so that subfunctions will be found
rmpath(genpath('.\_gsdata_'));

inputs.tail.orientation = 'centerline';  % either 'centerline' for original choice of aligned with body centerline at rear pole, or 'pole2pole' to align with a line passing through both polar semicircle centers

%load sweep_body2.mat %contains swept AR1 and AR2

inputs.bugtype = 'bacteria';  %bacteria with rotating tail, or dino with deforming tranverse and longitudinal flagella
inputs.do_timestepping = false;

switch inputs.bugtype
    case 'bacteria'
%         load(['../sweeps/',sweep_tempfile]);  %contains all new AR1, AR2, amp, lambda, nlambda
new_sweep.AR1 = 1;  new_sweep.AR2 = 0;  new_sweep.amp = NaN; new_sweep.lambda = NaN;  new_sweep.nlambda = NaN;
    case {'dino', 'sheet'}
        inputs.accuracy.check_intersections_tolerance = 0.15 / 2;  %hard code tail radius for now
        inputs.tail.motorBC = 'none';  %placeholder to appease mex compiler
end

%location to put output dumps in
%inputs.paths.dumpfolder = '../swept_dumps/';

%suffix to append to output filenames
inputs.paths.suffix = '';
%suffix = ['rotatetailangle_',num2str(rotate_tail_angle*180/pi),'_fast'];

%freeswim (set motorBC and have self-propelled motion)
% or
%forced (forced translation in x, y, z and/or forced rotation around x, y, z axes)
switch inputs.bugtype
    case 'bacteria'
        %         inputs.problemtype = {'forced','freeswim'}';
        inputs.problemtype = 'forced';
%         inputs.BC_type = [1];
    case 'dino'
        inputs.problemtype = {'freeswim'}';
%          inputs.potatohead = {[0 1 0 0], [0 1 0 1], [0 1 1 1], [0 0 0 1],   [0 0 1 0] , [0 1 1 0], [0 0 1 1],[1 1 0 0], [1 0 1 0],[1 1 0 1],[1 1 1 0],[1 1 1 1]}';  %potato head cases to do:  [Body Transverse Tail] on or off
%                  inputs.potatohead = {[1 0 1 1]}';  %potato head cases to do:  [Body Transverse Tail Wingtip] on or off
%      inputs.potatohead = {[1 0 1 0 0 0]}';  %potato head cases to do:  [Body Transverse Tail Coplanar_Hairs Normal_Top_Hairs Normal_Bottom_Hairs] on or off
%       inputs.potatohead = {[1 1 0 0 0 0],[1 0 0 1 0 0],[1 0 0 0 1 0],[1 0 0 0 0 1]}';
   inputs.potatohead = {[1 1 0 1 1 1],[1 1 0 1 0 0]}';
  inputs.potatohead = {[0 0 0 1 0 0] };
%   inputs.potatohead = {[1 1 1 1 1 1] };
  
%   inputs.potatohead = {[1 1 1 1 0 0] };

% inputs.BC_type = [1 1 1 1 1 1]; % 1 for no slip, 2 for free slip
    case 'sheet'
        inputs.problemtype = {'freeswim'};
        
end

%ignore hydrodynamic interaction integrals between body and tail(s)?
inputs.ignore_interaction = false;
switch getenv('computername')
    case 'UBERTOP'
        inputs.performance.nthreads = 8; %feature('numCores');  %number of threads for parallelization
    case {'CFD01','CFD02','CFD03','CFD04'}
        inputs.performance.nthreads = 20; %feature('numCores');  %number of threads for parallelization
end

inputs.performance.kinematics_interpolation = true;  %Use true as long as symmetry dictates that all kinematics in the body frame are periodic.  Must use false for problems with multiple bugs or walls that break symmetry.
inputs.performance.debug_mode = false;  %if true, *a lot* of additional output on adaptive surface integrals is saved and code is probably slower and definitely needs *a lot* of RAM
inputs.performance.randomize_verts = true;  %randomize order of verts within each submesh to reduce load balancing problem due to some body-body, tail-tail integrals taking forever
inputs.performance.numels_max = double(intmax('int32'));  %the real Coder limit may be smaller than the stated intmax?  % was 1E9 but could be due to a code error now fixed
inputs.performance.verbose = true;  %display timing info and other messages?

inputs.performance.timestepping.initial_length = 10000;
inputs.performance.timestepping.chunk_size = 10000;

%going for < 0.01 % error WRT all accuracy parameters

switch inputs.bugtype
    case 'bacteria'
        factor = 1;
    case {'dino', 'sheet'}
        factor = 2  ;
end

f2 = 1/2/2/2/2/2 ; % seems good
inputs.accuracy.epsilon =  4E-5 *factor ;  %just how regular are the regularized Stokeslets?  check convergence as epsilon --> 0

%all abs tolerances are later scaled by individual element surface area, so
%that big elements are allowed to have more total error than small
%elements, i.e.
%             these abs tolerances are per unit surface area!
% Note that since viscosity (mu) is accounted for in the solution after the
% matrix solve, it doesn't affect any integration tolerances, i.e. a
% viscosity of unity is effectively used during matrix assembly and solving

%traction tolerances are for integrating the regularized Stokeslet velocity
%field function (S) over the surface when forming the boundary integral
%equations
inputs.accuracy.integration_tol.traction.abstol =   16E-3   *factor  * f2 ; %scale by min (worst case) element size compared to convergence test case
inputs.accuracy.integration_tol.traction.reltol = 0;  %control error via abstol
inputs.accuracy.integration_tol.traction.maxevals = Inf; %don't limit maxevals

%force tolerances are for integrating simply 1 * hS * phi over the
%surface.  this is used primarily for integrating force, e.g. for the
%free-swimming force balance equations.
inputs.accuracy.integration_tol.force.abstol = 1E-4   *factor  * f2; %1E-7
inputs.accuracy.integration_tol.force.reltol = 0;
inputs.accuracy.integration_tol.force.maxevals = Inf;

%torque tolerances are for integrating r * hS * phi over the surface
%where r is a lever-arm from a specified reference point.  this is
%primarily used for integrating torques, e.g. for the free-swimming
%torque balance equations
inputs.accuracy.integration_tol.torque.abstol = 1E-4  *factor  * f2; %1E-5
inputs.accuracy.integration_tol.torque.reltol = 0;  %2.5E-5
inputs.accuracy.integration_tol.torque.maxevals = Inf;

switch inputs.bugtype
    case 'bacteria'
        inputs.accuracy.integration_tol.area.abstol = 1E-12 ;
        inputs.accuracy.integration_tol.area.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.integration_tol.area.maxevals = Inf;
    case {'dino', 'sheet'}
        inputs.accuracy.integration_tol.area.abstol = 1E-6 ;
        inputs.accuracy.integration_tol.area.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.integration_tol.area.maxevals = Inf;
end



%volume tolerances are for calculating mesh volume by a surface
%integral of F dot n via divergence theorem (F = [1/3 x, 1/3 y, 1/3 z])
switch inputs.bugtype
    case 'bacteria'
        inputs.accuracy.integration_tol.volume.abstol = 1E-12 ;
        inputs.accuracy.integration_tol.volume.reltol = 1E-1  ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.integration_tol.volume.maxevals = Inf;
    case 'dino'
        inputs.accuracy.integration_tol.volume.abstol = 1E-6 ;
        inputs.accuracy.integration_tol.volume.reltol = 1E-1  ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.integration_tol.volume.maxevals = Inf;
    case 'sheet'
        inputs.accuracy.integration_tol.volume.abstol = Inf;  %2D sheet, volume meaningless
        inputs.accuracy.integration_tol.volume.reltol = 1E-1  ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.integration_tol.volume.maxevals = Inf;
end



%centroid tolerances are for calculating centroid of enclosed volume
%integral of F dot n via divergence theorem (F = [1/2 x^2, 0, 0] or [0, 1/2 y^2, 0] or [0, 0, 1/2 z^2])
switch inputs.bugtype
    case 'bacteria'
        inputs.accuracy.integration_tol.centroid.abstol = 1E-12 ;
        inputs.accuracy.integration_tol.centroid.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.integration_tol.centroid.maxevals = Inf;
    case 'dino'
        inputs.accuracy.integration_tol.centroid.abstol = 1E-6 ;
        inputs.accuracy.integration_tol.centroid.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.integration_tol.centroid.maxevals = Inf;
    case 'sheet'
        inputs.accuracy.integration_tol.centroid.abstol = Inf;
        inputs.accuracy.integration_tol.centroid.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.integration_tol.centroid.maxevals = Inf;
end





% inputs.accuracy.interpolation.max_rel_diff = 0.1 ; %max allowable relative difference between last two most accurate trigonometric interpolants, over all interpolation fields (ordinarily U(1:3), Omega(1:3), omega)

% inputs.accuracy.interpolation.max_rel_diff = 0.5 ;
% inputs.accuracy.interpolation.max_rel_diff = 1 ;
% inputs.accuracy.interpolation.max_rel_diff = 0.5 ;

% inputs.accuracy.interpolation.max_rel_diff = 0.1 ;  % was 0.5 with old RMS fun diff

% inputs.accuracy.interpolation.max_rel_diff.mag = 0.005 ; % 0.05 often OK
% inputs.accuracy.interpolation.max_rel_diff.mag = 0.02 ;
% inputs.accuracy.interpolation.max_rel_diff.mag = 0.075 ;
% inputs.accuracy.interpolation.max_rel_diff.mag = 0.035 ;
% inputs.accuracy.interpolation.max_rel_diff.mag = 0.025 ;
inputs.accuracy.interpolation.max_rel_diff.mag = 0.1 / 100 ;
inputs.accuracy.interpolation.max_rel_diff.dir = 0.03 ;

inputs.accuracy.interpolation.vector_normalization = true;
inputs.accuracy.interpolation.max_iter = 5;  %max # adaptive interpolation iterations before giving up
% 1   2    3     4    5    6     7           max_iter
% 2   4    8     16   32   64    128          # phase angle evals


%criteria for determining convergence of avg swimming speed:
%(ode45 timestepping stops when this is reached)
% inputs.accuracy.timestepping.datafrac = 1; %last fraction of avg speed data to use in determining convergence(allows one to ignore initial big transient)
% inputs.accuracy.timestepping.cvtol = 1E-8;  %max coefficient of variation in last fraction of avg speed data allowable (5E-5 ~ 0.005% error expected)
% inputs.accuracy.timestepping.minpts = 50;  %min number of timesteps allowable (needed to prevent immediate fake convergence)
% inputs.accuracy.timestepping.error_tol = 0.0001;  %0.01 % error, where
% error is estimated as ( max(speeds) - min(speeds) ) / mean(speeds) / 2
% where speeds is all the avg speed data from the last interrogation interval time window.
% Factor of 2 is there since we can be very confident that the true avg
% speed is close to the avg value of the max and min speeds, since the
% trend seems to always be a decaying periodic signal.  Therefore this
% is still most likely a pretty conservative error estimate.
% inputs.accuracy.timestepping.normalized_interrogation_interval = 2.5E5 ;  %normalized to motor freq:  divide by freq to get interrogation interval
% ode45 error tolerances:
inputs.accuracy.timestepping.initialstep = 1E-6;
inputs.accuracy.timestepping.reltol = 1E-6;
inputs.accuracy.timestepping.abstol = 1E-12 ;
switch inputs.bugtype
    case 'bacteria'
inputs.accuracy.timestepping.T_interrogate = 2.^[-2:11];
inputs.accuracy.timestepping.diff_tols = [0.1:0.1:5];  %first entry is preferred max % diff between 2 consecutive avg speed estimates.  if this fails after T_interrogate(end), then tol is incremented and convergence checked again

    case 'dino'
       inputs.accuracy.timestepping.T_interrogate = 2.^[4:11];  % need enough time to generate decent swimming videos
       inputs.accuracy.timestepping.diff_tols = [0.05 0.1:0.1:5];  %first entry is preferred max % diff between 2 consecutive avg speed estimates.  if this fails after T_interrogate(end), then tol is incremented and convergence checked again

end
% inputs.accuracy.timestepping.diff_tols = [0.1:0.1:5];  %first entry is preferred max % diff between 2 consecutive avg speed estimates.  if this fails after T_interrogate(end), then tol is incremented and convergence checked again

switch inputs.bugtype
    case 'bacteria'
        inputs.accuracy.check_intersections_n_angles = 100;  %how many equally spaced phase angles to test for the body hitting the tail.  10 seems safe.  Nope, need more than 10:  trying 20....  Nope, 20 not always enough.  50?...
    case {'dino', 'sheet'}
        inputs.accuracy.check_intersections_n_angles = []; %not rotating tail like in bacteria case
end

inputs.constants.mu = 1E-3 / 1E6 * 1E9; %viscosity in kg / (s micron)   now in ug / (s micron), so forces will be in ug um s^-2 = pN
%gets multiplied by boundary integral
inputs.constants.multfactor = 1/8/pi;  %remove mu here and add back later to improve numerical scaling of A matrix
inputs.constants.power = 1E-3;  %power assumed to be dissipated by all bugs when calculating adjusted swimming speed

inputs.output.interpolation.doplot = true;  %plot of interpolation convergence progress
inputs.output.interpolation.doprint = true;  %save pdfs of current progres plot, if doplot == true
inputs.output.interpolation.fignum = 100;

inputs.output.timestepping.doplot = true;  %plot of timestepping convergence progress (for avg swimming speed)
inputs.output.timestepping.plotfreq = 10E3;  %update plot every plotfreq-th timestep
inputs.output.timestepping.doprint = true;  %save pdfs of current progress plot, if doplot == true
inputs.output.timestepping.fignum = 200 ;

if strcmp(inputs.bugtype, 'bacteria')
    
    %% Body Geometry
    
    % ellipsoid, capsule, curved_rod, dino
    %     inputs.body.shape = 'capsule';
    %     inputs.body.shape = 'ellipsoid';
    
    inputs.body.shape = 'curved_rod';
    
    
    inputs.body.suffix = '';
    inputs.body.V = 1;  %volume, microns^3
    % %equivalent sphere radius
    inputs.body.sphererad = (inputs.body.V*3/4/pi)^(1/3);
    
    
    
    % ARs = AR1(204);
    % AR2s = AR2(204);
    
    % ARs = 5.123;
    % AR2s = 0.543;  %to compare timing of orig vs memory saving version of matrix_assembly_mex
    
    % ARs = 2;  %seemingly largest body mesh
    % AR2s = 0.55; %seemingly largest tail mesh
    
    
    %node 195   1:20       21:40       41:60
    %node 196   61:80      81:100      101:120
    %node 197   121:140    141:160     161:180
    %node 198   181:194    195:208     209:224
    
    % ARs = AR1(61:120); %read from sweep_body.mat file at top of this script
    % AR2s = AR2(61:120);
    
    % ARs = AR1(209:224); %read from sweep_body.mat file at top of this script
    % AR2s = AR2(209:224);
    
    
    ARs = new_sweep.AR1;
    AR2s = new_sweep.AR2;
    
    if length(ARs) ~= length(AR2s)
        error('ARs and AR2s are corresponding vectors and must be same length.');
    end
    
    inputs.body.AR = mat2cell([ARs(:)'; AR2s(:)'],2,ones(1,length(ARs)));
    %each cell of ARcell is a different specific geometry to combine with other
    %swept parameters
    
    
    %% Tail
    
    %body only or also tail(s)?
    inputs.include_tail = false;
    if inputs.include_tail
        inputs.tail.tail_type = 'normal';  %either 'normal' for helical tail or 'debug' for dumb body-shaped "tail" with fewer elements
        if strcmp(inputs.tail.tail_type,'debug')  %override tail geometry parameters with a specific mesh file
            inputs.paths.debug_tail = 'ellipsoid_AR1_4_AR2_4';
        end
        %initially rotate tail by some angle to test for effects on friction coeffs, diffusivity in forced-flow simulations?
        inputs.tail.rotate_tail = false;
        inputs.tail.rotate_tail_angle = pi;
        
        % inputs.tail.suffix = {'_base'};
        inputs.tail.suffix = '';
        % Tail motor parameters
        
        %freq (set constant tail rotation rate, torque varies)
        %or
        %torque (set constant motor torque, rotation rate varies)
        % inputs.tail.motorBC = {'torque','freq'};
        inputs.tail.motorBC = 'torque';
        %in both cases, power theoretically varies
        inputs.tail.motor_torque = 1e-06;  %micron N microns   yields about 100 Hz rotation rate = 617 rad / sec
        % inputs.tail.motor_freq = 100 *2*pi ;  %rad / sec
        inputs.tail.motor_freq = 468 ;  %rad / sec
        
        
        
        % Tail Geometry
        switch inputs.tail.tail_type
            case 'normal'
                inputs.tail.radius = 0.05 * inputs.body.sphererad;  %ala Shum et al
                %                 inputs.tail.radius = 0.031018;
            case 'debug'
                inputs.tail.radius = 1.5633/2 + 0.1;
        end
        
        inputs.accuracy.check_intersections_tolerance = inputs.tail.radius / 2;
        
        %   factors = [0.75 1 1.25];
        % factors = 0.75;
        
        %helix wavelength, nondimensionalized as in Shum et al and then
        %redimensionalized...
        inputs.tail.lambda = num2cell(4.68     * inputs.body.sphererad * [1.25  1.5]);
        %   inputs.tail.lambda = 3.6291;
        %helix amplitude, nondimensionalized as in Shum et al and then
        %redimensionalized...
        %inputs.tail.amp = 0.87  / (2*pi/inputs.tail.lambda) ;
        inputs.tail.amp = num2cell(0.402 * [1.25   1.5]) ;
        % inputs.tail.amp = .5025;
        inputs.tail.nlambda = num2cell(1.49 * [ 0.875  1  1.125 ]);  %number of helix wavelengths in tail
        % inputs.tail.nlambda = 1.8625;
        %
        %    inputs.tail.lambda = new_sweep.lambda;
        %    inputs.tail.amp = new_sweep.amp;
        %    inputs.tail.nlambda = new_sweep.nlambda;
        
        inputs.tail.lambda = mat2cell([new_sweep.lambda(:)'],1,ones(1,length(new_sweep.lambda)));
        inputs.tail.amp = mat2cell([new_sweep.amp(:)'],1,ones(1,length(new_sweep.amp)));
        inputs.tail.nlambda = mat2cell([new_sweep.nlambda(:)'],1,ones(1,length(new_sweep.nlambda)));
        
    else
         inputs.tail.radius = NaN;  inputs.tail.amp = NaN;  inputs.tail.lambda = NaN;  inputs.tail.nlambda = NaN;    inputs.tail.suffix = '';
         inputs.tail.tail_type = 'normal';
        inputs.tail.motorBC = 'none';  %placeholder to appease mex compiler
    end
    
end
%% Parameter sweep

[fieldpaths, values] = traverse_inputs_struct(inputs); %convert nested struct to list of individual fieldpaths and values of each parameter (possibly cell arrays or vectors of values)

% separate out the cell arrays to be swept
swept_parameters = {};  swept_fieldpaths = {};
for i = 1:length(values)
    if iscell(values{i}) && size(values{i},2) == 1 %column vector cell array - should sweep over values of this parameter type
        swept_parameters{end+1} = values{i};
        swept_fieldpaths{end+1} = fieldpaths{i};
    end
end

% separate out the cell arrays to be iterated over one by one
iter_parameters = {};  iter_fieldpaths = {};
for i = 1:length(values)
    if iscell(values{i}) && size(values{i},1) == 1 %row vector cell array - should iterate over each matched group, without sweeping all possible combinations
        iter_parameters{end+1} = values{i};
        iter_fieldpaths{end+1} = fieldpaths{i};
    end
end

% form nested structure "paths" from lists of field names
swept_str = {};
for j = 1:length(swept_fieldpaths)
    swept_str{j} = [];
    for k = 1:length(swept_fieldpaths{j})
        swept_str{j} = [swept_str{j},'.',swept_fieldpaths{j}{k}];
    end
end

% form nested structure "paths" from lists of field names
clear iter_str
for j = 1:length(iter_fieldpaths)
    iter_str{j} = [];
    for k = 1:length(iter_fieldpaths{j})
        iter_str{j} = [iter_str{j},'.',iter_fieldpaths{j}{k}];
    end
end

sweep = parameter_sweep(swept_parameters);  %computes all combinations of each swept parameter type
%each row of sweep is a sweep combination, column is parameter type

%reorder sweep to alternate forced and freeswim, so that most of A matrix
%can be reused for each pair
if isequal(inputs.problemtype,{'forced','freeswim'}') || isequal(inputs.problemtype,{'freeswim','forced'}')
    % num rows of sweep should be factor of 2 since we have forced and
    % freeswim
    sweep2 = [];
    halflength = size(sweep,1) / 2;
    for i = 1:halflength
        sweep2 = [sweep2; sweep(i,:); sweep(i+halflength,:)];
    end
    sweep = sweep2;  clear sweep2
end

clear iter_parameters2
for i = 1:length(iter_parameters)
    temp = repmat(iter_parameters{i},size(sweep,1),1);
    iter_parameters2{i} = {temp{:}}';
end

iter_parameters2 = horzcat(iter_parameters2{:});  %iterated parameters, all runs

sweep2 = repmat(sweep,length(iter_parameters{1}),1);  %swept parameters, all runs

sweep = [sweep2  iter_parameters2];  %combine them for all iterated/swept parameters for all runs
str = [swept_str  iter_str];

%finally, make final Inputs struct, which has many indexed entries, each entry containing
%all the parameter values for a single job in the parameter sweep
clear Inputs
for i = 1:size(sweep,1) %each sweep iteration
    Inputs(i) = inputs; %lazily just copy all of inputs first, then overwrite swept fields
    for j = 1:size(sweep,2)  %swept parameter types
        eval(['Inputs(i)',str{j},' = sweep{i,j};']);
    end
end


% precompute a few heavily used things
for i = 1:length(Inputs)
    Inputs(i).accuracy.eps2 = Inputs(i).accuracy.epsilon.^2;
end



for i = 1:length(Inputs)
    switch inputs.bugtype
        case 'bacteria'
            %construct base name of output files
            Inputs(i).paths.namebase.body = [Inputs(i).body.shape,'_AR1_',num2str(Inputs(i).body.AR(1)),'_AR2_',num2str(Inputs(i).body.AR(2)),Inputs(i).body.suffix];
            switch Inputs(i).tail.tail_type
                case 'normal'
                    Inputs(i).paths.namebase.tail = ['tail_','radius_',num2str(Inputs(i).tail.radius,digits),'_amp_',num2str(Inputs(i).tail.amp,digits),'_lambda_',num2str(Inputs(i).tail.lambda,digits),'_nlambda_',num2str(Inputs(i).tail.nlambda,digits),Inputs(i).tail.suffix];
                case 'debug'
                    Inputs(i).paths.namebase.tail = [Inputs(i).paths.debug_tail];
            end
            
            %  rest = ['_eps_',num2str(Inputs(i).accuracy.epsilon),'_abstol_',num2str(Inputs(i).accuracy.integration_tol.traction.abstol),'_interptol_',num2str(Inputs(i).accuracy.interpolation.max_rel_diff),Inputs(i).paths.suffix];
            rest = Inputs(i).paths.suffix;
            
            switch Inputs(i).problemtype
                case 'freeswim'  %output file names use combined info of body and tail geometries
                    Inputs(i).paths.namebase.full = [Inputs(i).paths.namebase.body,'_',Inputs(i).paths.namebase.tail,'_motorBC_',Inputs(i).tail.motorBC, rest];
                case 'forced'
                    if Inputs(i).include_tail
                        Inputs(i).paths.namebase.full = [Inputs(i).paths.namebase.body,'_',Inputs(i).paths.namebase.tail,'_forced',rest];
                    else
                        Inputs(i).paths.namebase.full = [Inputs(i).paths.namebase.body,'_forced',rest];
                    end
            end
            
        case 'dino'
  
            suf = [];
            names = {'Body','Transverse','Tail','Coplanar_Hairs','Normal_Top_Hairs','Normal_Bottom_Hairs'};
            last_present = find(Inputs(i).potatohead == 1,1,'last');
            for n = 1:length(names)
                if Inputs(i).potatohead(n)
                    if n ~= last_present
                        suf = [suf , names{n} , '-' ];
                    else
                        suf = [suf , names{n} ];
                    end
                end
            end
            
            %   prefix = temp_sweep{i};
            if ~isempty(prefix)
                Inputs(i).paths.namebase.full = [prefix,'_',suf];
            else
                Inputs(i).paths.namebase.full = [suf];
            end
            
            Inputs(i).paths.interp_partial_file = [Inputs(i).paths.namebase.full , '.interp_partial.mat'];
            
            
        case 'sheet'
            Inputs(i).paths.namebase.full = 'sheet3';
    end
    
end




%create stuff if it isn't there
% if ~exist(inputs.paths.intersections_file,'file')
%     intersections = [];
%     save(inputs.paths.intersections_file,'intersections');
% end

if ~exist(inputs.paths.dumpfolder,'dir')
    mkdir(inputs.paths.dumpfolder);
end

if ~exist(inputs.paths.results_folder,'dir')
    mkdir(inputs.paths.results_folder);
end

if ~exist([inputs.paths.results_folder,inputs.paths.results_file],'file')
    Results = [];
    save([inputs.paths.results_folder,inputs.paths.results_file],'Results');
end

clear sweep* iter_parameters* inputs values fieldpaths
