%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear geom

body_radius = 15.2;
body_height0 = 40;  % was 10.7
scale_x = 0.5;  % was 0.825

geom.body.center = [body_height0 * scale_x / 2 + body_height0 / 2    0      0]';

geom.phase_speed = 2*pi*46;  %if lowest common period T = 4*pi s, covered in 2*pi rad, is 1/2 rad/sec for the entire beat cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geom.tail.radius = 0.15   * 1;
% geom.tail.radius = 0.9;


geom.tail.lambda = 22.2;
geom.tail.nlambda = 1.8531 ;
geom.tail.amp = [1.752*0.95  11.4 / 2];
% geom.tail.kE = [3 * 0.24      2 * 0.05 * 1.5 ];
geom.tail.kE = [3 * 0.24  * 0.5     2 * 0.05 * 1.5 ];


geom.tail.omega = 2 * pi * 46 ;   %rad/sec
% period = 2*pi rad / (2*pi*46  rad/sec) = 1/46 sec
% phase speed = 2*pi / period = 2*pi*46  rad/sec

% geom.tail.t_transition = 3.050645083561573;
% geom.tail.t_transition = 2.999558294702207;
% geom.tail.t_transition = 3.000076881165709;
geom.tail.t_transition = 4.173965824914811 ; % with offset
geom.tail.t_transition = 3.319438443750901;  % no offset

% geom.tail.translation = [8.326846940519667 0.000000000000004 14.100000000000000]';
% geom.tail.translation = [8.146345068104413 0.000000000000004 14.650000000000000]';
% geom.tail.translation = [30.250000000000028 0.000000000000004 14.650000000000000]';
% geom.tail.translation = [25.182751825774218 0.000000000000004 14.650000000000000]';
geom.tail.translation = [30.019250099175309 0.000000000000004 14.650000000000000]';  % with offset
geom.tail.translation = [27.000000000000004 0.000000000000004 14.650000000000000]'; % no offset

geom.tail.t_min = 0;
geom.tail.t_max = 2*pi * geom.tail.nlambda;  %t value where hemispherical end cap is centered
safety_factor = 0.1;  % t value at start of startign sph or end of ending sphere must be within (t_max - t_min)*safety_factor of t_min or t_max respectively

%geom.tail_angle = -25 *0  ;

geom.tail.rotation_angle = [-pi/2 ]';  %if multiple values, several consecutive rotations are done
% -90 deg vs + 90 deg in Salome because in Salome order of rotation is
% x then z, while here the z rotation is hardcoded into tail eqs and
% occurs first
geom.tail.rotation_pt  = [0 0 0]';
geom.tail.rotation_vec = [1 0 0]';

% if we are doing centered_tail:
geom.tail.shift = [11 0 0]' - geom.tail.translation;  % x = 11 seems to work for fat tail, need to check that 2nd part matches current Parameters file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geom.transverse.R = 13.45;  %orig before wingtip
% geom.transverse.R = 13.45 + 0.5;  %shallow groove for wingtip
geom.transverse.R = body_radius + 9.2 - 10;

%geom.R = 15.7;  big groove I think
geom.transverse.d = 0;
%geom.w = 0.1;
geom.transverse.w = 0.85;  %4
geom.transverse.rf = 3.2 / 2;
geom.transverse.c = 18;
geom.transverse.omega = 2 * pi * 46;

% geom.transverse.b = 3.7 ;  %wavelength of underlying helix - 0 for no helix, just a circle
% geom.transverse.b = 5 ;  %wavelength of underlying helix - 0 for no
% helix, just a circle     with offset
geom.transverse.b = 0 ;  % no offset
geom.transverse.N_revs = 1  ;   % * 0.3
%     geom.u_max = 2*pi*geom.N_revs;
% geom.transverse.u_max = 6.130414606345004; %no longer based on N_revs since sheet stops where groove ends, which Salome determines
% geom.transverse.u_max = 6.126525207868735;
geom.transverse.u_max = 6.130084811492812  ;  % with offset
geom.transverse.u_max = 6.283185307179586 ;


geom.transverse.v_max = 1;
%     geom.shift = 8.7;
% geom.transverse.shift = body_height0 * 0.825 * 1.25;
geom.transverse.shift = 32.5;  %gotten from Salome command line, need to add shift_x to Parameters.txt     with offset
geom.transverse.shift = 30.000000000000004;  % no offset
%%%%%%%%%%%%%%%%%%%%%%
% geom.transverse.wingtip_type = 'parallel_2';
geom.hairs.Coplanar_Hairs.h_min = 0;  % always 0 for parallel wingtip
geom.hairs.Coplanar_Hairs.h_max = 2;

geom.hairs.Normal_Top_Hairs.h_min = -1.5;  
geom.hairs.Normal_Top_Hairs.h_max = 0;  %top is actually negative side

geom.hairs.Normal_Bottom_Hairs.h_min = 0;  
geom.hairs.Normal_Bottom_Hairs.h_max = 1.5;  %bottom is actually positive side

% geom.transverse.h_min = -2;  %appears to be upper edge for normal wingtip
% geom.transverse.h_max = 2;  %appears to be lower edge for normal wingtip
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%