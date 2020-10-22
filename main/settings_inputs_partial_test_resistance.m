% shat 2 get working as mex:  save_Col_inds_mexed , store mesh constants,
% matrix assembly mexed, compute_BCs_mexed

% also, investigate 2% difference in field_vel from
% field_velocity.m vs field_velocity_mexed

% improve T6interp by using simplified formulas if all alpha, beta, gamma
% for the Mesh = 0.5 (precompute this once beforehand and send in a flag).

% check if shat is converged vs integration tol for field velocity for dino

% check if shat converged vs basic integration tols for dino !! make sure
% to test at least 2 consecutive refinements since original choices seem to
% be a rel min of error in traction for a forced sphere

% test error due to neglecting double layer potential by comparing known
% velocity at surfaces with velocity computed there via velocity field
% calc?  as per Dave around ~eq 12 in nearest neighbor pape, and cilia
% section of Dave RSBEM pape ?

% change epsilon to be sheet thickness for transverse and hair sheets?

% might be possible to speedup by evaluating some or all of the 54 SL
% integrals per element per collocation pt separately, but would have to
% avoid redoing same calcs already done in calcS?  Unknown how much # evals
% varies among the 54 integrals.

% IMPORTANT
 %%%%%%%%%%%% At some point need to figure out constants in front of all terms e.g. 1/8/pi
 
 
% IMPORTANT
% check speed of mexed version with / without if/else shortcuts in calcS,
% calcT - one would expect shortcut to be faster, but calcT without shortcut is
% ~5 times faster unmexed

% check various speed optimizations in integrand, calcS, calcT, e.g. is it still better to just compute all entries even if we know things are symmetric?

% recycling of A for bacteria resistance / mobility problem combo will likely be broken
% now due in part to rejiggering of velocity BCs

% shortcut for bacterial mobility where A matrix entries are rotated at
% each timestep is prolly buggered due to adding DL 

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
        inputs.paths.datfolder = 'C:\Hull\shum test\';  %refined version in \sphere mesh\ %location of .dat T6 mesh file(s)
        inputs.paths.datfolder = 'C:\Hull\test meshes\';
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
inputs.kinematics_interp_method = "spline";  % either "trig" for Fourier interpolation (won't work if any interp phases are discarded!) or "spline" (still works with discarded points)
inputs.periodic_kinematics = true;
assert(all(~ismember(roundn([0 2*pi],inputs.interp_phases_tol),roundn(inputs.discard_interp_phases,inputs.interp_phases_tol))) , 'Not allowed to discard phase 0 or 2*pi since this is where periodicity is enforced');

addpath(genpath('./')); %add all subfolders to path so that subfunctions will be found
rmpath(genpath('.\_gsdata_'));

inputs.Tail.orientation = "centerline";  % either "centerline" for original choice of aligned with body centerline at rear pole, or "pole2pole" to align with a line passing through both polar semicircle centers

%load sweep_body2.mat %contains swept AR1 and AR2

inputs.bugtype = "bacteria";  %bacteria with rotating tail, or dino with deforming tranverse and longitudinal flagella
inputs.rotating_flagellum = false;  % does the problem involve a rotating tail, i.e. bacteria?  
inputs.do_timestepping = false;

switch inputs.bugtype
    case "bacteria"
        %         load(['../sweeps/',sweep_tempfile]);  %contains all new AR1, AR2, amp, lambda, nlambda
        new_sweep.AR1 = 1;  new_sweep.AR2 = 0;  new_sweep.amp = 0.402; new_sweep.lambda = 2.9032;  new_sweep.nlambda = 1.49;
    case {"dino", "sheet"}
        inputs.accuracy.check_intersections_tolerance = 0.15 / 2;  %hard code tail radius for now
        inputs.Tail.motorBC = "none";  %placeholder to appease mex compiler
end

%location to put output dumps in
%inputs.paths.dumpfolder = '../swept_dumps/';

%suffix to append to output filenames
inputs.paths.suffix = '';
%suffix = ['rotatetailangle_',num2str(rotate_tail_angle*180/pi),'_fast'];

%mobility (set motorBC and have self-propelled motion)
% or
%resistance (forced translation in x, y, z and/or forced rotation around x, y, z axes)

inputs.BC_type.default = 1; % 1 for no slip, 2 for free slip

inputs.performance.eliminate_DL.default = false; % eliminate DL integral from BIE:  should always do this for sheets, for which it is not even clear how to integrate the DL
% and it always cancels out
% for rigid bodies, don't have to eliminate the DL, but should to save CPU time, since there is no loss of accuracy.
% for deforming bodies, should not eliminate the DL.  If it is eliminated, one is solving for the traction jump across the boundary.  For resistance
% problems, this means that there is no way to compute drag force (which should just be the integral of the outer traction), though velocity fields should be right.  
% For mobility problems, a force-free constraint will still work and swimming trajectories and velocity fields should be right, but one can't compute drag
% on individual bodies in multiple-body systems (e.g. body + tail), nor can one compute power dissipation over the boundaries.
% this is under the performance field even though it could also belong under accuracy if true for a deforming / frees slip body...

switch inputs.bugtype
    case "bacteria"
        %         inputs.problemtype = {"resistance","mobility"}';
        inputs.problemtype = "mobility";
        inputs.problemtype = "resistance";
%         inputs.BC_type = [1];

inputs.BC_type.Body = 1;
inputs.BC_type.Tail = 1;

inputs.parent_topology.Body = "closed";  % topology of overall parent body that this submesh is a part of, e.g. if a closed surface is split into multiple 
% submeshes, parent_toplogy = "closed" even though individual submeshes appear to be "open" on their own
inputs.parent_topology.Tail = "closed"; % topology is used when calculating solid angle - for an open surface, there is "outside" fluid on both sides and 
% thus solid angle is 4 pi; for intersections of a closed and open surface, the open surface doesn't affect solid angle at the intersection

inputs.performance.eliminate_DL.Body = true;
inputs.performance.eliminate_DL.Tail = true;

inputs.performance.DL_singularity_removal.Body = true; 
inputs.performance.DL_singularity_removal.Tail = true;

inputs.is_mesh_rigid.Body = true; % is mesh rigid (only translates/rotates) over time?  Note that body may be rigid but if it is remeshed over time, mesh isn't rigid
inputs.is_mesh_rigid.Tail = true; % boundary integrals over rigid submeshes can be quickly computed for mobility problems by rotating initial tensor values instead of redoing integrals at every phase / time

% below, outer cell needed just to tell code not to try to do a parameter sweep over this
% inputs.coincident_submeshes = { {"Body" , "Tail"  } };  % each cell contains a physically separate group of coincident submeshes that share nodes
% if no submeshes are coincident, each submesh should be a separate cell here
% if some closed geometric entities are split into submeshes (e.g. to facilitate calculations of force on different regions), this is designated using a 2nd
% level of cells.  e.g. {  {["Body1", "Body2"], "Transverse"} , "Appendage", { ["Tail1", "Tail2"] }  } indicates 3 non-coincident entities:  the first
% consists of two submeshes forming one closed surface, the body, and the transverse sheet sharing nodes with the body.  the second is some other appendage
% completely separate in space.  the third is another completely separate closed surface, the tail, consisting of two submeshes.  Specifying this info
% allows us to correctly deal with a) creation of global node indices across geometrically isolated, independently meshed entities and b) calculation of solid angles for which all submeshes comprising a closed
% surface must be handled as a unit
% an "entity" is a single closed or open surface, e.g. here entities are Body1+Body2, Transverse, Appendage, Tail1+Tail2.  
inputs.coincident_submeshes = { { {"Body" } , {"Tail"} } };  %need gratuitous outer {} here just to avoid trying to do a parameter sweep
inputs.coincident_submeshes = { { {"Body" } } }; 


    case "dino"
        inputs.problemtype = {"mobility"}';
%          inputs.potatohead = {[0 1 0 0], [0 1 0 1], [0 1 1 1], [0 0 0 1],   [0 0 1 0] , [0 1 1 0], [0 0 1 1],[1 1 0 0], [1 0 1 0],[1 1 0 1],[1 1 1 0],[1 1 1 1]}';  %potato head cases to do:  [Body Transverse Tail] on or off
%                  inputs.potatohead = {[1 0 1 1]}';  %potato head cases to do:  [Body Transverse Tail Wingtip] on or off
%      inputs.potatohead = {[1 0 1 0 0 0]}';  %potato head cases to do:  [Body Transverse Tail Coplanar_Hairs Normal_Top_Hairs Normal_Bottom_Hairs] on or off
%       inputs.potatohead = {[1 1 0 0 0 0],[1 0 0 1 0 0],[1 0 0 0 1 0],[1 0 0 0 0 1]}';
   inputs.potatohead = {[1 1 0 1 1 1],[1 1 0 1 0 0]}';
  inputs.potatohead = {[0 0 0 1 0 0] };
%   inputs.potatohead = {[1 1 1 1 1 1] };
  
%   inputs.potatohead = {[1 1 1 1 0 0] };

inputs.BC_type.Body = 1;
inputs.BC_type.Transverse = 1;
inputs.BC_type.Tail = 1;
inputs.BC_type.Coplanar_Hairs = 1;
inputs.BC_type.Normal_Top_Hairs = 1;
inputs.BC_type.Normal_Bottom_Hairs = 1;

inputs.parent_topology.Body = "closed";
inputs.parent_topology.Transverse = "open";
inputs.parent_topology.Tail = "closed";
inputs.parent_topology.Coplanar_Hairs = "open";
inputs.parent_topology.Normal_Top_Hairs = "open";
inputs.parent_topology.Normal_Bottom_Hairs = "open";

inputs.performance.eliminate_DL.Body = false;
inputs.performance.eliminate_DL.Transverse = true;
inputs.performance.eliminate_DL.Tail = false;
inputs.performance.eliminate_DL.Coplanar_Hairs = true; % the DL always cancels out for sheets, there one always solves for the traction jump across the sheet
inputs.performance.eliminate_DL.Normal_Top_Hairs = true;
inputs.performance.eliminate_DL.Normal_Bottom_Hairs = true;

% whether to use singularity removal identity "trick" (Poz Practical Guide to BEM p. 179) when evaluating DL integrals that we're not eliminating entirely
% the identity is only valid for closed surfaces but this should always be the case?  (for sheets, the DL integral is always eliminated or zero so this switch doesn't matter)
% typically one would completely eliminate the DL for rigid bodies so this switch is mainly for deforming bodies, where DL should be kept if we care about
% power, and we should use the trick to speed up code
inputs.performance.DL_singularity_removal.Body = true; 
inputs.performance.DL_singularity_removal.Transverse = NaN;
inputs.performance.DL_singularity_removal.Tail = true;
inputs.performance.DL_singularity_removal.Coplanar_Hairs = NaN;
inputs.performance.DL_singularity_removal.Normal_Top_Hairs = NaN;
inputs.performance.DL_singularity_removal.Normal_Bottom_Hairs = NaN;




inputs.is_mesh_rigid.Body = false;  % you'd think true, but currently it is remeshed every flagellar phase due to being coincident with transverse, so the mesh isn't rigid
inputs.is_mesh_rigid.Transverse = false;
inputs.is_mesh_rigid.Tail = false;
inputs.is_mesh_rigid.Coplanar_Hairs = false;
inputs.is_mesh_rigid.Normal_Top_Hairs = false;
inputs.is_mesh_rigid.Normal_Bottom_Hairs = false;


inputs.coincident_submeshes = { { "Body", "Transverse", "Coplanar_Hairs", "Normal_Top_Hairs", "Normal_Bottom_Hairs" }, {"Tail"} };
% each uppermost level cell contains a physically separate group of coincident submeshes that share nodes
% if no submeshes are coincident, each submesh should be in a separate cell
% within each upper level cell, the inner cells contain geometrically separate entities.  The main point of this is to specify closed surfaces that might
% be split into multiple (open) submeshes for post-processing convenience (e.g. to faciliate calculations of force on different regions)
% e.g. { { ["Body1", "Body2"], "Transverse", "Coplanar_Hairs", "Normal_Top_Hairs", "Normal_Bottom_Hairs" }, {"Tail"} } 
% indicates 3 upper-level non-coincident groups:  the first consists of the body, and a few open surfaces that all share nodes.  the body is further split into 
% two submeshes.  the second upper-level group contains just the tail, completely disconnected in space from the body/sheets group.
% Specifying this info allows us to correctly deal with a) creation of global node indices across geometrically isolated, independently meshed entities and 
% b) calculation of solid angles for which all submeshes comprising a single closed surface must be handled as a unit
% since sheets are assumed to always have solid angle = 4 pi, it doesn't matter whether sheets made of multiple coincident submeshes are considered as
% single geometric entities or not - the distinction is only important for closed surfaces

    case "sheet"
        inputs.problemtype = {"mobility"};
        
end


%ignore hydrodynamic interaction integrals between body and tail(s)?
% Dave suggests the matrix built with this true might be sparse and good for preconditioning the full system (with this false) with an iterative solver
inputs.accuracy.mesh.ignore_interaction = false;

inputs.performance.rigid_body_matrix_rotation = false;  % if there are any rigid bodies in a free-swimming problem, e.g. cell body or rigid bacterial tail, there is no 
% need to recompute boundary integral contributions at new times or beat phases for collocation pts and elements both on the same rigid body.  Instead, can apply a rotation to the
% initial tensors / vectors in the matrix / RHS for each BIE in the system.  Can also do this with combo resistance/mobility runs on the same geometry, 
% computing everything for the resistance problem and then recycling the rigid body integrals for all times/phases in the mobility problem.
% This should be faster than recomputing everything, though it does require more memory to store the recycled initial values in addition to the current
% full matrix.  This setting is used in conjunction with the "rigid" flag of each submesh - note that the *mesh* itself must be rigid over time, so we 
% can't use this shortcut for rigid bodies that are nonetheless remeshed over time (e.g. dino body, currently).
% This switch may later get set to false in assembly_input if we don't actually have any rigid submeshes in the problem
inputs.performance.rigid_body_tensor_storage = "sparse";  % if above switch is true, this controls whether the recycled values are stored as a "sparse" vs
% "full" array.  "sparse" should often be better but can try testing.

switch getenv('computername')
    case 'UBERTOP'
        inputs.performance.nthreads = 8; %feature('numCores');  %number of threads for parallelization
    case {'CFD01','CFD02','CFD03','CFD04'}
        inputs.performance.nthreads = 20; %feature('numCores');  %number of threads for parallelization
end

inputs.performance.kinematics_interpolation = true;  %Use true as long as symmetry dictates that all kinematics in the body frame are periodic.  Must use false for problems with multiple bugs or walls that break symmetry.
inputs.performance.debug_mode = false;  %if true, *a lot* of additional output on adaptive surface integrals is saved and code is probably slower and definitely needs *a lot* of RAM
inputs.performance.randomize_nodes = false;  %randomize order of verts within each submesh to reduce load balancing problem due to some body-body, tail-tail integrals taking forever
% instead of randomizing meshes, try randomizing just the parfor assembly loop (see comment in matrix_assembly_mex?)
inputs.performance.numels_max = intmax('int32');  %the real Coder limit may be smaller than the stated intmax?  % was 1E9 but could be due to a code error now fixed
inputs.performance.verbose = true;  %display timing info and other messages?

inputs.performance.timestepping.initial_length = 10000;
inputs.performance.timestepping.chunk_size = 10000;

%going for < 0.01 % error WRT all accuracy parameters

%just how regular are the regularized Stokeslets?  check convergence as epsilon --> 0

switch inputs.bugtype
    case "bacteria"
        factor = 1;
    case {"dino", "sheet"}
        factor = 2  ;
end

% f2 = 1/2/2/2/2/2 ;
f2 = 1/2/2/2    *100  ;

inputs.accuracy.mesh.epsilon.default =  4E-5 *factor; % used if a value not specified for a particular submesh
inputs.accuracy.mesh.epsilon.Body = 4E-5; inputs.accuracy.mesh.epsilon.Tail = 4E-5;
inputs.accuracy.mesh.epsilon.Transverse = 10E-9 * 1E6 / 2;  % transverse sheet is theoretically 2 layers of cell membrane, each 5 nm thick.  eps should be the radius of a thin structure
inputs.accuracy.mesh.epsilon.Coplanar_Hairs = 10E-9 * 1E6 / 2; % conveniently, hairs all also theoretically 10 nm wide
inputs.accuracy.mesh.epsilon.Normal_Top_Hairs = 10E-9 * 1E6 / 2;
inputs.accuracy.mesh.epsilon.Normal_Bottom_Hairs = 10E-9 * 1E6 / 2;

inputs.accuracy.network.epsilon = 1/4/pi;
% it is unclear what happens for verts shared between submeshes, e.g. Body
% and Transverse - would need to review code that decides unique global
% verts, but prolly doesn't matter much either way


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
inputs.accuracy.mesh.integration_tol.stokeslet.abstol =   16E-3   *factor  * f2 ; %scale by min (worst case) element size compared to convergence test case
inputs.accuracy.mesh.integration_tol.stokeslet.reltol = 0;  %control error via abstol
inputs.accuracy.mesh.integration_tol.stokeslet.maxevals = Inf; %don't limit maxevals

%force tolerances are for integrating simply 1 * hS * phi over the
%surface.  this is used primarily for integrating force, e.g. for the
%free-swimming force balance equations.
inputs.accuracy.mesh.integration_tol.force.abstol = 1E-4   *factor  * f2; %1E-7
inputs.accuracy.mesh.integration_tol.force.reltol = 0;
inputs.accuracy.mesh.integration_tol.force.maxevals = Inf;

%torque tolerances are for integrating r * hS * phi over the surface
%where r is a lever-arm from a specified reference point.  this is
%primarily used for integrating torques, e.g. for the free-swimming
%torque balance equations
inputs.accuracy.mesh.integration_tol.torque.abstol = 1E-4  *factor  * f2; %1E-5
inputs.accuracy.mesh.integration_tol.torque.reltol = 0;  %2.5E-5
inputs.accuracy.mesh.integration_tol.torque.maxevals = Inf;

switch inputs.bugtype
    case "bacteria"
        inputs.accuracy.mesh.integration_tol.area.abstol = 1E-12 ;
        inputs.accuracy.mesh.integration_tol.area.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.mesh.integration_tol.area.maxevals = Inf;
    case {"dino", "sheet"}
        inputs.accuracy.mesh.integration_tol.area.abstol = 1E-6 ;
        inputs.accuracy.mesh.integration_tol.area.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.mesh.integration_tol.area.maxevals = Inf;
end



%volume tolerances are for calculating mesh volume by a surface
%integral of F dot n via divergence theorem (F = [1/3 x, 1/3 y, 1/3 z])
switch inputs.bugtype
    case "bacteria"
        inputs.accuracy.mesh.integration_tol.volume.abstol = 1E-12 ;
        inputs.accuracy.mesh.integration_tol.volume.reltol = 1E-1  ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.mesh.integration_tol.volume.maxevals = Inf;
    case "dino"
        inputs.accuracy.mesh.integration_tol.volume.abstol = 1E-6 ;
        inputs.accuracy.mesh.integration_tol.volume.reltol = 1E-1  ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.mesh.integration_tol.volume.maxevals = Inf;
    case "sheet"
        inputs.accuracy.mesh.integration_tol.volume.abstol = Inf;  %2D sheet, volume meaningless
        inputs.accuracy.mesh.integration_tol.volume.reltol = 1E-1  ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.mesh.integration_tol.volume.maxevals = Inf;
end



%centroid tolerances are for calculating centroid of enclosed volume
%integral of F dot n via divergence theorem (F = [1/2 x^2, 0, 0] or [0, 1/2 y^2, 0] or [0, 0, 1/2 z^2])
switch inputs.bugtype
    case "bacteria"
        inputs.accuracy.mesh.integration_tol.centroid.abstol = 1E-12 ;
        inputs.accuracy.mesh.integration_tol.centroid.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.mesh.integration_tol.centroid.maxevals = Inf;
    case "dino"
        inputs.accuracy.mesh.integration_tol.centroid.abstol = 1E-6 ;
        inputs.accuracy.mesh.integration_tol.centroid.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.mesh.integration_tol.centroid.maxevals = Inf;
    case "sheet"
        inputs.accuracy.mesh.integration_tol.centroid.abstol = Inf;
        inputs.accuracy.mesh.integration_tol.centroid.reltol = 1E-1 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        inputs.accuracy.mesh.integration_tol.centroid.maxevals = Inf;
end

% double layer regularized stresslet integrals

% integrals of plain reg. stresslet T.n, i.e. for U.T.n term
% inputs.accuracy.integration_tol.stresslet.Tn.abstol = 1E-6;
% inputs.accuracy.integration_tol.stresslet.Tn.reltol = 0;
% inputs.accuracy.integration_tol.stresslet.Tn.maxevals = Inf;

% integrals of r * T.n, i.e. for (Omega x r).T.n term
% inputs.accuracy.integration_tol.stresslet.rTn.abstol = 20 * inputs.accuracy.integration_tol.stresslet.Tn.abstol; % typical value of r for dino is prolly 20 um
% inputs.accuracy.integration_tol.stresslet.rTn.reltol = 0;
% inputs.accuracy.integration_tol.stresslet.rTn.maxevals = Inf;

% integrals of phi * T.n, i.e. for u.T.n term
% inputs.accuracy.integration_tol.stresslet.phiTn.abstol = 1E-6; 
inputs.accuracy.mesh.integration_tol.stresslet.phiTn.abstol = 1E-2; 
inputs.accuracy.mesh.integration_tol.stresslet.phiTn.reltol = 0;
inputs.accuracy.mesh.integration_tol.stresslet.phiTn.maxevals = Inf;



% % integrals of v * reg. stresslet . n, i.e. for v.T.n term where v is known
% % velocity (of e.g. flagella) relative to the body moving at U + Omega x r
% % or, in resistance problems where v is known velocity of body in fixed frame
% inputs.accuracy.integration_tol.stresslet.plain.abstol = 500* inputs.accuracy.integration_tol.stresslet.plain.abstol; % edge of transverse sheet moves at around 500 um/s
% inputs.accuracy.integration_tol.stresslet.plain.reltol = 0;
% inputs.accuracy.integration_tol.stresslet.plain.maxevals = Inf;


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
    case "bacteria"
inputs.accuracy.timestepping.T_interrogate = 2.^[-2:11];
inputs.accuracy.timestepping.diff_tols = [0.1:0.1:5];  %first entry is preferred max % diff between 2 consecutive avg speed estimates.  if this fails after T_interrogate(end), then tol is incremented and convergence checked again

    case "dino"
       inputs.accuracy.timestepping.T_interrogate = 2.^[4:11];  % need enough time to generate decent swimming videos
       inputs.accuracy.timestepping.diff_tols = [0.05 0.1:0.1:5];  %first entry is preferred max % diff between 2 consecutive avg speed estimates.  if this fails after T_interrogate(end), then tol is incremented and convergence checked again

end
% inputs.accuracy.timestepping.diff_tols = [0.1:0.1:5];  %first entry is preferred max % diff between 2 consecutive avg speed estimates.  if this fails after T_interrogate(end), then tol is incremented and convergence checked again

switch inputs.bugtype
    case "bacteria"
        inputs.accuracy.check_intersections_n_angles = 100;  %how many equally spaced phase angles to test for the body hitting the tail.  10 seems safe.  Nope, need more than 10:  trying 20....  Nope, 20 not always enough.  50?...
    case {"dino", "sheet"}
        inputs.accuracy.check_intersections_n_angles = []; %not rotating tail like in bacteria case
end

inputs.constants.mu = 1E-3 / 1E6 * 1E9; %viscosity in kg / (s micron)   now in ug / (s micron), so forces will be in ug um s^-2 = pN
%gets multiplied by boundary integral
inputs.constants.multfactor = 1/8/pi;  %remove mu here and add back later to improve numerical scaling of A matrix
inputs.constants.power = 1E-3;  %power assumed to be dissipated by all bugs when calculating adjusted swimming speed


inputs.accuracy.triangle_integration.reference_nodes = [0 1 0; 0 0 1];  %reference triangle normalized node coords
inputs.accuracy.triangle_integration.rule = 2;
[inputs.accuracy.triangle_integration.G, inputs.accuracy.triangle_integration.Weights, inputs.accuracy.triangle_integration.PTS] =  SMPRMS( 2, inputs.accuracy.triangle_integration.rule );  %always 2 dimensions and constant rule # for all integrals (2nd arg is rule #)





inputs.output.interpolation.doplot = true;  %plot of interpolation convergence progress
inputs.output.interpolation.doprint = true;  %save pdfs of current progres plot, if doplot == true
inputs.output.interpolation.fignum = 100;

inputs.output.timestepping.doplot = true;  %plot of timestepping convergence progress (for avg swimming speed)
inputs.output.timestepping.plotfreq = 10E3;  %update plot every plotfreq-th timestep
inputs.output.timestepping.doprint = true;  %save pdfs of current progress plot, if doplot == true
inputs.output.timestepping.fignum = 200 ;

if strcmp(inputs.bugtype, "bacteria")
    
    %% Body Geometry
    
    % ellipsoid, capsule, curved_rod, dino
    %     inputs.body.shape = "capsule";
%         inputs.Body.shape = "ellipsoid";
    
    inputs.Body.shape = "curved_rod";
    
    
    inputs.Body.suffix = '';
    inputs.Body.V = 1;  %volume, microns^3
    % %equivalent sphere radius
    inputs.Body.sphererad = (inputs.Body.V*3/4/pi)^(1/3);
 
    ARs = new_sweep.AR1;
    AR2s = new_sweep.AR2;
    
    if length(ARs) ~= length(AR2s)
        error('ARs and AR2s are corresponding vectors and must be same length.');
    end
    
    inputs.Body.AR = mat2cell([ARs(:)'; AR2s(:)'],2,ones(1,length(ARs)));
    %each cell of ARcell is a different specific geometry to combine with other
    %swept parameters
    
    
    %% Tail
    
    %body only or also tail(s)?
    inputs.include_tail = false;
    if inputs.include_tail
        inputs.Tail.tail_type = "normal";  %either "normal" for helical tail or "debug" for dumb body-shaped "tail" with fewer elements
        if strcmp(inputs.Tail.tail_type,"debug")  %override tail geometry parameters with a specific mesh file
            inputs.paths.debug_tail = 'ellipsoid_AR1_4_AR2_4';
        end
        %initially rotate tail by some angle to test for effects on friction coeffs, diffusivity in resistance problem simulations?
        inputs.Tail.rotate_tail = false;
        inputs.Tail.rotate_tail_angle = pi;
        
        % inputs.tail.suffix = {'_base'};
        inputs.Tail.suffix = '';
        % Tail motor parameters
        
        %freq (set constant tail rotation rate, torque varies)
        %or
        %torque (set constant motor torque, rotation rate varies)
        % inputs.tail.motorBC = {"torque","freq"};
        inputs.Tail.motorBC = "torque";
        %in both cases, power theoretically varies
        inputs.Tail.motor_torque = 1E3;  % ug um/s^2 um = fN um (femto N um)  yields about 100 Hz rotation rate = 617 rad / sec
        % max E coli torque reported to be 2000 pN nm = 2000 fN um (KK Mandadapu - ‎2015) so this works out
        
        % inputs.tail.motor_freq = 100 *2*pi ;  %rad / sec
        inputs.Tail.motor_freq = 468 ;  %rad / sec
        
        
        
        % Tail Geometry
        switch inputs.Tail.tail_type
            case "normal"
%                 inputs.Tail.radius = 0.05 * inputs.Body.sphererad;  %ala Shum et al
                                inputs.Tail.radius = 0.031018;
                                inputs.Tail.radius = 0.08;
            case "debug"
                inputs.Tail.radius = 1.5633/2 + 0.1;
        end
        
        inputs.accuracy.check_intersections_tolerance = inputs.Tail.radius / 2;
        
        %   factors = [0.75 1 1.25];
        % factors = 0.75;
        
        %helix wavelength, nondimensionalized as in Shum et al and then
        %redimensionalized...
        inputs.Tail.lambda = num2cell(4.68     * inputs.Body.sphererad * [1.25  1.5]);
        %   inputs.tail.lambda = 3.6291;
        %helix amplitude, nondimensionalized as in Shum et al and then
        %redimensionalized...
        %inputs.tail.amp = 0.87  / (2*pi/inputs.tail.lambda) ;
        inputs.Tail.amp = num2cell(0.402 * [1.25   1.5]) ;
        % inputs.Tail.amp = .5025;
        inputs.Tail.nlambda = num2cell(1.49 * [ 0.875  1  1.125 ]);  %number of helix wavelengths in tail
        % inputs.tail.nlambda = 1.8625;
        %
        %    inputs.tail.lambda = new_sweep.lambda;
        %    inputs.tail.amp = new_sweep.amp;
        %    inputs.tail.nlambda = new_sweep.nlambda;
        
        inputs.Tail.lambda = mat2cell([new_sweep.lambda(:)'],1,ones(1,length(new_sweep.lambda)));
        inputs.Tail.amp = mat2cell([new_sweep.amp(:)'],1,ones(1,length(new_sweep.amp)));
        inputs.Tail.nlambda = mat2cell([new_sweep.nlambda(:)'],1,ones(1,length(new_sweep.nlambda)));
        
    else
         inputs.Tail.radius = NaN;  inputs.Tail.amp = NaN;  inputs.Tail.lambda = NaN;  inputs.Tail.nlambda = NaN;    inputs.Tail.suffix = '';
         inputs.Tail.tail_type = "normal";
        inputs.Tail.motorBC = "none";  %placeholder to appease mex compiler
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

%reorder sweep to alternate resistance and mobility problems, so that most of A matrix
%can be reused for each pair
if isequal(inputs.problemtype,{"resistance","mobility"}') || isequal(inputs.problemtype,{"mobility","resistance"}')
    % num rows of sweep should be factor of 2 since we have resistance and
    % mobility
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


for i = 1:length(Inputs)
    switch inputs.bugtype
        case "bacteria"
            %construct base name of output files
            Inputs(i).paths.namebase.Body = Inputs(i).Body.shape + '_AR1_' + num2str(Inputs(i).Body.AR(1)) + '_AR2_' + num2str(Inputs(i).Body.AR(2)) + Inputs(i).Body.suffix;
            switch Inputs(i).Tail.tail_type
                case "normal"
                    Inputs(i).paths.namebase.Tail = "Tail_" + "radius_" + num2str(Inputs(i).Tail.radius,digits) + "_amp_" + num2str(Inputs(i).Tail.amp,digits) + "_lambda_" + num2str(Inputs(i).Tail.lambda,digits) + "_nlambda_" + num2str(Inputs(i).Tail.nlambda,digits) + Inputs(i).Tail.suffix;
                case "debug"
                    Inputs(i).paths.namebase.Tail = [Inputs(i).paths.debug_tail];
            end
            
            %  rest = ['_eps_',num2str(Inputs(i).accuracy.epsilon),'_abstol_',num2str(Inputs(i).accuracy.integration_tol.traction.abstol),'_interptol_',num2str(Inputs(i).accuracy.interpolation.max_rel_diff),Inputs(i).paths.suffix];
            rest = Inputs(i).paths.suffix;
            
            switch Inputs(i).problemtype
                case "mobility"  %output file names use combined info of Body and tail geometries
                    Inputs(i).paths.namebase.full = Inputs(i).paths.namebase.Body + '_' + Inputs(i).paths.namebase.Tail + '_motorBC_' + Inputs(i).Tail.motorBC +  rest;
                case "resistance"
                    if Inputs(i).include_tail
                        Inputs(i).paths.namebase.full = Inputs(i).paths.namebase.Body + '_' + Inputs(i).paths.namebase.Tail + '_resistance' + rest;
                    else
                        Inputs(i).paths.namebase.full = Inputs(i).paths.namebase.Body + '_resistance' + rest;
                    end
            end
            
        case "dino"
  
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
            
            
        case "sheet"
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


% Inputs.paths.namebase.Body = "curved_rod_AR1_1_AR2_0"
% Inputs.paths.namebase.Tail = "short_tail"
% Inputs.paths.namebase.full = "sphere_short_tail"
% Inputs.Body.shape = "ellipsoid"

% Inputs.paths.namebase.Body = "ellipsoid_AR1_1.67_AR2_1.67"
Inputs.paths.namebase.Body = "curved_rod_AR1_1_AR2_0"
% Inputs.paths.namebase.Tail = "tail_radius_0.031018_amp_0.402_lambda_2.9032_nlambda_1.49"
Inputs.paths.namebase.Tail = "longer_tail"
Inputs.paths.namebase.Tail = "short_tail2"
Inputs.paths.namebase.full = "shum_test"
% Inputs.Body.shape = "ellipsoid"
