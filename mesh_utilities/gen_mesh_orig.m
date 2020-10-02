%generates body and tail CAD geometry and meshes via Salome
clear geom mesh paths
digits = 16;
%geom.shape = 'ellipsoid'; %triaxial ellipsoid
%geom.shape = 'capsule';   %elliptic capsule (cross section can be an ellipse)
%   geom.shape = 'curved_rod'; %capsule bent around a circle
geom.shape = 'tail';   %helical flagellum capped with spheres
%geom.shape = 'dino';  %idealized dinoflagellate body, cylindrical center section and ellipsoidal ends
init_only = false;
%init_only = true;
% redo_failures_only = true;
redo_failures_only = false;

paths.outfolder = 'E:/Hull/swept_meshes/';  %be sure to use / instead of \ here; regexprep won't work with \
paths.infolder = 'E:/Hull/salome/';

outputprefix = geom.shape;  %put this at beginning of output file names
suffix = '';  %append at end of file names, if used, start with _

pyname = [paths.infolder,outputprefix,'_sweep',suffix,'_',num2str(randi(intmax))];  %python script that goes into Salome with random number for uniqueness

%OscarKatie E.coli is ~ 1.9 um long and 0.9 um wide, giving volume ~ 1 um^3
V = 1; %all body shapes will have this volume = 1 um^3
sphererad = (V*3/4/pi)^(1/3);  %equivalent sphere radius

if ~strcmp(geom.shape,'tail')
    geom.V = V ;  %all body shapes will have this volume = 1 um^3
    geom.sphererad = sphererad;  %equivalent sphere radius
end

%% compute geometry dimensions based on input parameters and aspect ratios

switch geom.shape
    case 'tail'
        
        %points_per_turn = [200]; %how many discretized points along splined centerline per helix wavelength
        % factors = [0.75 1 1.25];
        %factors = [1.25];
        % factors = 1;
        %number of wavelengths for entire tail
        nlambda = 1.49 * [ 0.875  1  1.125 ];  %power optimized, Shum et al
        %  geom.nlambda = 1.28;  %torque optimized, Shum et al
        
        %wavelength of helix
        lambda = [4.68 * [1.25  1.5]  ] * sphererad; %power opt
        %geom.lambda = 1.4 * geom.sphererad; %torque opt
        
        %amplitude of helix i.e. radius of helix centerline
        %  amp = [0.87  ./ (2*pi/lambda)];   %power opt
        amp = [0.402  * [1.25   1.5]  ];   %power opt
        %amp = 0.68 / (2*pi/geom.lambda);  %torque opt
        
        %  pipeRadius = [0.05] * sphererad;  %radius of flagellum "pipe"
        % in Fluid Mech of Propulsion by cilia and flagella p. 352, range
        % of flagella radius stated to be 0.012 -  0.02 um.  Shum used
        % about 0.03 = 0.05 * sphererad
        pipeRadius_orig = [0.05] * sphererad;
        
        
        
        [Nlambda, Lambda, Amp] = ndgrid(nlambda, lambda, amp);
        
        Nlambda = Nlambda(:);  Lambda = Lambda(:);  Amp = Amp(:);
        
        
        
        %         Nlambda = new_tails(start:last,3);  Lambda = new_tails(start:last,2);  Amp = new_tails(start:last,1);
        %
        %           Nlambda = new_tails(:,3);  Lambda = new_tails(:,2);  Amp = new_tails(:,1);
        
        start = 1;  %66 fuxored
        last = length(new_sweep.AR1);
        last = 1;
        
        Nlambda = new_sweep.nlambda(start:last)';  Lambda = new_sweep.lambda(start:last)';  Amp = new_sweep.amp(start:last)';
        
        
        pipeRadius = repmat(pipeRadius_orig,length(Nlambda),1); %avg value of stated range
        
        % kE = k in Shum et al
        kE = 2 * pi ./ Lambda;
        
        
        for c = 1:length(Nlambda)
            geom(c) = geom(1); %copy defaults / constants
            paths(c) = paths(1);
            
            geom(c).nlambda = Nlambda(c);
            geom(c).lambda = Lambda(c);
            geom(c).amp = Amp(c);
            geom(c).pipeRadius = pipeRadius(c);
            geom(c).kE = kE(c);
            
            %wavenumber = k = 2 * pi / lambda
            % t = k * xi
            % xi = lambda * t / 2 / pi
            % 1/k = lambda / 2 / pi
            % k = 2 * pi / lambda
            factor = 1 *0.9; %1.375
            %mesh(c).maxsize = 0.08 * geom(c).pipeRadius / pipeRadius_orig * factor;
            mesh(c).maxsize = 0.065 * geom(c).pipeRadius / pipeRadius_orig * factor;
            % mesh(c).maxsize = 0.056;
            mesh(c).minsize = mesh(c).maxsize / 2;
            mesh(c).fineness = 4;  %from 0 to 4 for coarse to fine, and 5 for custom
            
            factor  = 1*0.9;
            mesh(c).maxsize_ends = 0.01875 * geom(c).pipeRadius / pipeRadius_orig * factor;  %0.015
            mesh(c).minsize_ends = mesh(c).maxsize_ends * 0.8;
            mesh(c).fineness_ends = 3; %from 0 to 4 for coarse to fine, and 5 for custom
            
        end
        
    case {'ellipsoid'}
        %3 radii a, b, c for ellipsoids (microns)
        
        geom.AR1 = 0.6; % a/b
        
        geom.AR2 = 7;   % a/c
        %aspect_ratio2 = ARtemp / aspect_ratio; % b/c
        
        
        geom.b = (3/4*geom.V / pi * geom.AR2 / geom.AR1^2 )^(1/3);
        geom.a = geom.b * geom.AR1;
        geom.c = geom.b * geom.AR1 / geom.AR2;
        
        
        
        mesh.maxsize = min([geom.a geom.b geom.c]) * 0.8   *0.4  *5 ;  % add /2 for mincurvfrac >= 0.3
        mesh.minsize = min([geom.a geom.b geom.c]) * 0.5   *0.4 *2 ; %0.9;  % /11
        
        mesh.fineness = 0;  %from 0 to 4 for coarse to fine, and 5 for custom
        %%
    case 'capsule'
        
        geom.AR1 = 0.6; % (total length) / (2 * major radius)  .4
        
        geom.AR2 = 7;  % (total length) / (2 * minor radius)  10
        
        %   aspect_ratio2 = ARtemp / aspect_ratio; %ellipsoids:  b/c
        
        %compute equivalent ellipsoid of rev
        %         b = (3/4*V / pi / aspect_ratio * aspect_ratio2 )^(1/3);
        %         a = b * aspect_ratio;
        %         c = b / aspect_ratio2;
        
        temp = roots([2/3*pi*(3*geom.AR1^2/geom.AR2 - (geom.AR2/geom.AR1)^(-4/3)) 0 0 -geom.V]); %see tablet worksheet
        for i = 1:length(temp)
            if isreal(temp(i)) && temp(i) > 0
                geom.major_radius = temp(i);  % major radius of cross section ellipse
                break
            end
        end
        
        if sum(imag(temp) == 0) > 1
            error('Problem with solving for geometry parameters');
        end
        
        geom.minor_radius = geom.major_radius / (geom.AR2/geom.AR1);  % minor radius of cross section ellipse
        
        geom.caps_c = (geom.major_radius^2 * geom.minor_radius)^(1/3);  % c minor radius of end cap ellipsoids
        
        geom.height = 2*geom.major_radius*(geom.AR1 - (geom.AR2/geom.AR1)^(-1/3));  %height of cylinder
        if geom.height < 0
            error(['Calculations yield negative height = ',num2str(geom.height)]);
        end
        geom.totlength = geom.height + 2*geom.caps_c;  %length of complete capsule
        
        mesh.maxsize = min([geom.major_radius,geom.minor_radius,geom.caps_c])  * .4 * 20;%0.62; %140 seems coarsest possible was 140 *.5
        mesh.minsize = min([geom.major_radius,geom.minor_radius,geom.caps_c])  * .16 *.9;%0.9;  %80 seems coarsest possible was 47   /6
        mesh.fineness = 0;  %from 0 to 4 for coarse to fine, and 5 for custom.  Recommend just using 0 - 4
        
        % maxsize = min(a,b)*0.7;%0.62; %140 seems coarsest possible was 140 *.5
        % minsize = min(a,b)*0.6;%0.9;  %80 seems coarsest possible was 47   /6
        %
        % fineness = 5;  %from 0 to 4 for coarse to fine, and 5 for custom in which case below parameters are used.  "coarse" seems to mean a large allowable range in element size
        % growth_rate = 0.10;  %fraction of neighboring element size an element is allowed to grow by  0.4
        % seg_per_edge = 12;  %doesn't seem to matter for ellipsoids
        % seg_per_radius = 1;  %1 seems to work alright, can be < 1 or > 1, < 1 tends to fail
        
        %         growth_rate = 0.3;  %fraction of neighboring element size an element is allowed to grow by  0.4
        %         seg_per_edge = 15;  %doesn't seem to matter for ellipsoids
        %         seg_per_radius = 2;  %1 seems to work alright, can be < 1 or > 1, < 1 tends to fail
        %
        %%
    case 'dino'
        
        geom.body.AR1 = 1.4;
        geom.body.AR2 = 0.5;
        
        geom.body.points_per_turn = 200; %how many discretized points along splined centerline per helix wavelength
        
        %number of wavelengths for entire tail
        geom.body.nlambda = 1.49;  %power optimized, Shum et al
        %  geom.nlambda = 1.28;  %torque optimized, Shum et al
        
        %wavelength of helix
        geom.lambda = 4.68 * geom.sphererad; %power opt
        %geom.lambda = 1.4 * geom.sphererad; %torque opt
        
        %amplitude of helix i.e. radius of helix centerline
        geom.amp = 0.87 / (2*pi/geom.lambda);   %power opt
        %amp = 0.68 / (2*pi/geom.lambda);  %torque opt
        
        geom.pipeRadius = 0.05 * geom.sphererad;  %radius of flagellum "pipe"
        
        geom.limit_tol = 0.00025  * 0;  %for fuse operation of spheres with helix, but maybe not anymore?...
        
        temp = roots([2*pi*(geom.AR1 - geom.AR2) + 4/3*pi*geom.AR2, 0, 0, -geom.V]);
        for i = 1:length(temp)
            if isreal(temp(i)) && temp(i) > 0
                geom.radius = temp(i);  % radius of circular cross section
                break
            end
        end
        
        geom.caps_c = geom.radius * geom.AR2;
        geom.height = 2*geom.radius*(geom.AR1 - geom.AR2);
        geom.totlength = geom.height + 2*geom.caps_c;
        
        mesh.maxsize = geom.radius * 0.7 * 0.7 / 2;
        mesh.minsize = geom.radius * 0.25 * 5 / 2 ;
        mesh.fineness = 0;
        
    case 'curved_rod'
        % AR1 = arclength / (2*radius)     arclength = total arclength
        % AR2 = arclength / (2*pi*radius of curv)
        
        AR1 = [1.05 1.1 1.5  2:(1):12];
        AR2 = 0.05:(0.05):0.95;
        
        
        AR1 = [2.5 3.5 4.5 5.5 6.5];
        AR2 = [0.425:0.05:1];
        
        
        [AR1, AR2] = ndgrid(AR1,AR2);
        AR1 = AR1(:);  AR2 = AR2(:);
        
        AR1 = [2.5      2.5   2.5   3.5   3.5    3.5   4.5   4.5   4.5   5.5   5.5   5.5   6.5    6.5     8.5  8.5  8.5  8.5  8.5  8.5  8.5 8.5 8.5 9.5 9.5 9.5 9.5 9.5 9.5 9.5 9.5 10.5 10.5  10.5 10.5 10.5 10.5 10.5 10.5 10.5 11.5  11.5 11.5 11.5 11.5 11.5 11.5 11.5]';
        AR2 = [0.35     0.15  0.05  0.35  0.25   0.05  0.35  0.15  0.05  0.35  0.15  0.05  0.25   0.15    0.1  0.2  0.3  0.4  0.5  0.6  0.7 0.8 0.9 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.1  0.2   0.3  0.4  0.5  0.6  0.7  0.8  0.9  0.1   0.3  0.4  0.5   0.6 0.7  0.8   0.9]';
        %         var1 = eff.curved_rod.AR1s;   var2 = eff.curved_rod.AR2s;
        %         X_Y = [[var1]; [var2] ]';  %two vars to use for X and Y axes
        %         [X_Y_unq,inds,inds2] = unique(X_Y,'rows');  %we already did these AR1, AR2 combos so don't redo meshes
        %
        %         AR1AR2 = [AR1 AR2];
        %         temp = setdiff(AR1AR2, X_Y_unq,'rows');  %remove combos we already have meshes for
        %         AR1 = temp(:,1);  AR2 = temp(:,2);
        AR1 = [1.75 1.75 6.5 10.5 8.5];
        AR2 = [0.025 0.125 0.075 0.85 0.85];
        
        
        AR1 = [2.7 2.9      4.75   4.75   4.75    5.25   5.25   5.25];
        AR2 = [0   0         0.625  0.65   0.675   0.625  0.65  0.675];
        
        %  AR1 = 1.01;  AR2 = 0.35;
        
        factor = 4.2 * 1;  %1 for base refinement
        
        invalids = false(size(AR1));  %will keep track of shapes that are actually impossible due to being too fat to form a donut
        
        for c = 1:length(AR1)
            geom(c) = geom(1); %copy defaults / constants
            paths(c) = paths(1);
            
            geom(c).AR1 = AR1(c);
            geom(c).AR2 = AR2(c);
            
            if geom(c).AR1 == 1  %we have a sphere
                geom(c).radius = geom.sphererad;
                geom(c).height = 0;
                geom(c).totlength = 2*geom.sphererad;
                geom(c).nturns = NaN;
                geom(c).radius_curv = NaN;
                
                
            elseif geom(c).AR2 == 0  % don't have a sphere, but do have a straight rod of zero curvature
                
                temp = roots([2/3*pi*(3*geom(c).AR1 - 1) 0 0 -geom(c).V]); %see tablet worksheet
                for i = 1:length(temp)
                    if isreal(temp(i)) && temp(i) > 0
                        geom(c).radius = temp(i);  % major radius of cross section ellipse
                        break
                    end
                end
                
                if sum(imag(temp) == 0) > 1
                    error('Problem with solving for geometry parameters');
                end
                
                geom(c).height = 2*geom(c).radius*(geom(c).AR1 - 1);  %height of cylinder
                if geom(c).height < 0
                    error(['Calculations yield negative height = ',num2str(geom(c).height)]);
                end
                geom(c).totlength = geom(c).height + 2*geom(c).radius;  %length of complete capsule
                
                %         mesh.maxsize = geom(c).radius  * .4 * 20;%0.62; %140 seems coarsest possible was 140 *.5
                %         mesh.minsize = geom(c).radius  * .16 *.9;%0.9;  %80 seems coarsest possible was 47   /6
                %         mesh.fineness = 0;  %from 0 to 4 for coarse to fine, and 5 for custom.  Recommend just using 0 - 4
                geom(c).radius_curv = Inf;  %this is how we inform Salome that this is a straight rod
                geom(c).nturns = NaN;
                
            else % have a general curved rod
                
                temp = roots([2*pi*(geom(c).AR1 - 1/3) 0 0 -geom(c).V]);
                for i = 1:length(temp)
                    if isreal(temp(i)) && temp(i) > 0
                        geom(c).radius = temp(i);
                        break
                    end
                end
                
                geom(c).arclength = 2*geom(c).radius*geom(c).AR1;
                geom(c).radius_curv = geom(c).arclength / geom(c).AR2 / 2 / pi;
                
                if geom(c).radius_curv - geom(c).radius/2 < eps  %if radius of inner boundary is zero or negative, this shape is impossible
                    invalids(c) = true;
                    disp('donut is too fat')
                    fateroo
                end
                
                closest_dist = 2*geom(c).radius_curv*sin( (geom(c).arclength - 2*geom(c).radius) / 2 / geom(c).radius_curv ) - 2*geom(c).radius; %closest distance between surfaces of end spheres
                if closest_dist < eps && geom(c).AR2 >= 0.5 % this check doesn't mean much for not uber curved rods because end spheres might just be close due to shortness of inner segment
                    % invalids(c) = true;
                    disp('end spheres might touch')
                    pause
                end
                
                geom(c).nturns = (geom(c).arclength - 2 * geom(c).radius) / 2 / pi / geom(c).radius_curv;  %will always be a fraction - using modified helix-creating code to generate curved rods in a plane
                geom(c).height = NaN;
                
                
            end
            
            
            %  factor = 0.3;  %0.155 gives just over the max limit of verts
            mesh(c).maxsize = ( geom(c).radius ^ (0.4) ) * 0.6  * factor    ;  %0.8
            mesh(c).minsize = mesh(c).maxsize * 1/2   *0.15          ;  %1/3

            mesh(c).fineness = 3;
            %             mesh(c).fineness = 1;
            
            %   mesh(c).fineness = 1;
            
            %too many total verts:  22612  11492  10952
            %works:  8688 9456  10628
            
            % 10600 verts works
            % 31800 equations works
            % 1.0112E9 matrix elements works
            
        end
        
        geom = geom(~invalids);  mesh = mesh(~invalids); paths = paths(~invalids);   %get rid of invalid shapes, don't even try meshing them
        
        
end


%% make filenames
for c = 1:length(geom)
    
    switch geom(c).shape
        
        case 'tail'
            
            namebase = [outputprefix,'_','radius','_',num2str(geom(c).pipeRadius,digits),'_','amp','_',num2str(geom(c).amp,digits),'_','lambda','_',num2str(geom(c).lambda,digits),...
                '_','nlambda','_',num2str(geom(c).nlambda,digits)];
            
        case {'ellipsoid', 'capsule','curved_rod','dino'}
            
            namebase = [outputprefix,'_','AR1','_',num2str(geom(c).AR1,digits),'_','AR2','_',num2str(geom(c).AR2,digits)];
            
    end
    
    
    meshname = [paths(c).outfolder,namebase,suffix]; %mesh dat file that comes out of Salome
    
    paths(c).mesh_file = [meshname,'.dat'];
    paths(c).metafile = [meshname,'_metadata','.mat'];  %metafile holds geometry and mesh parameters that go with dat file
    paths(c).python_file = [pyname,'.py'];
end

%%
if redo_failures_only && ~init_only  %only attempt to mesh cases that failed before.  keep existing meshes that worked.
    worked = false(1,length(geom));
    for c = 1:length(geom)
        c/length(geom)
        [~, Metadata] = load_mesh(paths(c).mesh_file, [], [], true);  %only load Metadata to check whether meshing failed before
        if ~isempty(Metadata) && isfield(Metadata.mesh,'meshing_succeeded') && Metadata.mesh.meshing_succeeded  %Metadata file exists and we checked the mesh and it appeared to have worked
            worked(c) = true;
        end
    end
    geom = geom(~worked);  mesh = mesh(~worked);  paths = paths(~worked);
end


%%

if ~init_only
    % delete existing mesh .dat files since Salome is too dumb to overwrite
    % them
    for c = 1:length(geom)
        if exist(paths(c).mesh_file,'file') %Salome 7.4.0 can't seem to overwrite an existing dat file
            delete(paths(c).mesh_file)
        end
    end
    
else  %if init_only, don't do anything else
    return
end

template = [paths(c).infolder,geom(c).shape,'.py'];  %each shape has template python script that we then insert parameter choices into
fid = fopen(template);
text = fread(fid, inf, '*char')';
fclose(fid);

%% modify template python scripts by inserting parameter values
fs = '%0.16g, ';  %formatting specification for writing values to text files
% insert commas between outputnames and get into necessary format
temp = [{paths.mesh_file}; repmat({', '},1,length(paths))];
pathsprint = [];
for c = 1:length(paths)
    pathsprint = [pathsprint, '''', paths(c).mesh_file, ''', '];
end

switch geom(1).shape  %all shapes in parameter sweep *should* be the same!
    case 'tail'
        if batch
            batch_input = 'True';
            done_filename = [pathsprint(1:end-6) 'done'''];
        else
            batch_input = 'False';
            done_filename = '''''';
        end
        text = regexprep(text, 'batch_input', batch_input  );
        % text = regexprep(text, 'points_per_turn', num2str(geom(c).points_per_turn));
        text = regexprep(text, 'amp_list_input',  ['[',sprintf(fs,[geom.amp]),']'] );
        text = regexprep(text, 'pipeRadius_list_input',  ['[',sprintf(fs,[geom.pipeRadius]),']'] );
        text = regexprep(text, 'wavelength_list_input',  ['[',sprintf(fs,[geom.lambda]),']'] );
        text = regexprep(text, 'nlambda_list_input',  ['[',sprintf(fs,[geom.nlambda]),']'] );
        text = regexprep(text, 'kE_list_input',  ['[',sprintf(fs,[geom.kE]),']'] );
        
        
        text = regexprep(text, 'maxsize_list_input',  ['[',sprintf(fs,[mesh.maxsize]),']'] );
        text = regexprep(text, 'minsize_list_input',  ['[',sprintf(fs,[mesh.minsize]),']'] );
        text = regexprep(text, 'fineness_list_input',  ['[',sprintf(fs,[mesh.fineness]),']'] );
        
        text = regexprep(text, 'maxsize_ends_list_input',  ['[',sprintf(fs,[mesh.maxsize_ends]),']'] );
        text = regexprep(text, 'minsize_ends_list_input',  ['[',sprintf(fs,[mesh.minsize_ends]),']'] );
        text = regexprep(text, 'fineness_ends_list_input',  ['[',sprintf(fs,[mesh.fineness_ends]),']'] );
        
        text = regexprep(text, 'outputname_list_input', pathsprint );
        text = regexprep(text, 'done_filename', done_filename );
        
        
        
    case 'ellipsoid'
        %make ellipsoid by starting with sphere and then scaling
        sphere_rad = min([geom(c).a geom(c).b geom(c).c]);  %not to be confused with equivalent sphere radius
        scale1 = geom(c).a/sphere_rad;
        scale2 = geom(c).b/sphere_rad;
        scale3 = geom(c).c/sphere_rad;
        
        text = regexprep(text, 'sphererad', num2str(sphere_rad,digits));
        text = regexprep(text, 'scale1, scale2, scale3', [num2str(scale1,digits),', ',num2str(scale2,digits),', ',num2str(scale3,digits)]);
        text = regexprep(text, 'outdat', paths(c).mesh_file);
        text = regexprep(text, 'maxsize', num2str(mesh(c).maxsize,digits));
        text = regexprep(text, 'minsize', num2str(mesh(c).minsize,digits));
        text = regexprep(text, 'fineness', num2str(mesh(c).fineness,digits));
        
        %         if fineness == 5
        %             text = regexprep(text, 'growthrate', num2str(growth_rate));
        %             text = regexprep(text, 'segperedge', num2str(seg_per_edge));
        %             text = regexprep(text, 'segperradius', num2str(seg_per_radius));
        %         end
        
    case 'capsule'
        % shape will start as cylindrical capsule, where radius = minor_radius
        scale1 = geom(c).caps_c / geom(c).minor_radius; %how much to stretch or compress in x direction to get c right
        height0 = geom(c).height / scale1;  %start height as this, then it becomes actual height after scaling
        scale2 = geom(c).major_radius / geom(c).minor_radius;
        scale3 = 1;  %keep minor_radius as z radius
        
        text = regexprep(text, 'Vertex_1, Vector_1, minor_radius, height0', ['Vertex_1, Vector_1, ',num2str(geom(c).minor_radius,digits),', ',num2str(height0,digits)] );
        text = regexprep(text, 'Cylinder_1, deltax, deltay, deltaz', ['Cylinder_1, -',num2str(height0/2,digits),', 0, 0']);
        text = regexprep(text, 'Vertex_2, minor_radius', ['Vertex_2, ',num2str(geom(c).minor_radius,digits)]);
        text = regexprep(text, 'Vertex_3, minor_radius', ['Vertex_3, ',num2str(geom(c).minor_radius,digits)]);
        text = regexprep(text, 'scale1, scale2, scale3', [num2str(scale1,digits),', ',num2str(scale2,digits),', ',num2str(scale3,digits)]);
        text = regexprep(text, 'outdat', paths(c).mesh_file);
        text = regexprep(text, 'maxsize', num2str(mesh(c).maxsize,digits));
        text = regexprep(text, 'minsize', num2str(mesh(c).minsize,digits));
        text = regexprep(text, 'fineness', num2str(mesh(c).fineness,digits));
        %text = regexprep(text, 'growthrate', num2str(growth_rate));
        %text = regexprep(text, 'segperedge', num2str(seg_per_edge));
        %text = regexprep(text, 'segperradius', num2str(seg_per_radius));
        
    case 'dino'
        % shape will start as cylindrical capsule
        scale1 = geom(c).caps_c / geom(c).radius; %how much to stretch or compress in x direction to get c right
        height0 = geom(c).height / scale1;  %start height as this, then it becomes actual height after scaling
        scale2 = 1; %keep cross section circular
        scale3 = 1;  %keep radius as z radius
        
        text = regexprep(text, 'Vertex_1, Vector_1, body_radius, body_height0', ['Vertex_1, Vector_1, ',num2str(geom(c).radius,digits),', ',num2str(height0,digits)] );
        text = regexprep(text, 'Cylinder_1, deltax, deltay, deltaz', ['Cylinder_1, -',num2str(height0/2,digits),', 0, 0']);
        text = regexprep(text, 'Vertex_2, radius', ['Vertex_2, ',num2str(geom(c).radius,digits)]);
        text = regexprep(text, 'Vertex_3, radius', ['Vertex_3, ',num2str(geom(c).radius,digits)]);
        text = regexprep(text, 'scale1, scale2, scale3', [num2str(scale1,digits),', ',num2str(scale2,digits),', ',num2str(scale3,digits)]);
        text = regexprep(text, 'outdat', paths(c).mesh_file);
        text = regexprep(text, 'maxsize', num2str(mesh(c).maxsize,digits));
        text = regexprep(text, 'minsize', num2str(mesh(c).minsize,digits));
        text = regexprep(text, 'fineness', num2str(mesh(c).fineness,digits));
        
    case 'curved_rod'
        
        text = regexprep(text, 'pipeRadius_list_input', ['[',sprintf(fs,[geom.radius]),']'] );
        text = regexprep(text, 'radius_list_input', ['[',sprintf(fs,[geom.radius_curv]),']'] );
        text = regexprep(text, 'nturns_list_input', ['[',sprintf(fs,[geom.nturns]),']'] );
        text = regexprep(text, 'height_list_input', ['[',sprintf(fs,[geom.height]),']'] );
        text = regexprep(text, 'maxsize_list_input', ['[',sprintf(fs,[mesh.maxsize]),']'] );
        text = regexprep(text, 'minsize_list_input', ['[',sprintf(fs,[mesh.minsize]),']'] );
        text = regexprep(text, 'fineness_list_input', ['[',sprintf(fs,[mesh.fineness]),']'] );
        
        text = regexprep(text, 'outputname_list_input', pathsprint );
        text = strrep(text,'NaN','float(''NaN'')');  text = strrep(text,'Inf','float(''Inf'')');
        
end



if exist(paths(1).python_file)
    delete(paths(1).python_file);
end
fid = fopen(paths(1).python_file, 'w');
fwrite(fid, text);
fclose(fid);



%% do the deed - run Salome and create mesh (Salome saves the mesh .dat file itself)
switch geom(1).shape
    case 'curved_rod'  %newest Salome version doesn't like curved rods.  sigh.
        [gar,bage] = system(['E:\Hull\salome\SALOME-7.7.1-WIN64\run_salome.bat -t ',paths(1).python_file]);
    case 'tail'
        [gar,bage] = system(['E:\Hull\salome\SALOME-7.4.0-WIN64\run_salome.bat -t ',paths(1).python_file]);  %was 7.4.0
        %         [gar,bage] = system(['E:\Hull\salome\SALOME-7.7.1-WIN64\run_salome.bat -t ',paths(1).python_file]);  %was 7.4.0
    otherwise
        [gar,bage] = system(['E:\Hull\salome\SALOME-7.5.1-WIN64\run_salome.bat -t ',paths(1).python_file]);
end
%% save metadata file that has geometry parameters that go with the mesh file
%meshname = dat_file;
for c = 1:length(geom)
    clear Metadata
    Metadata.geom = geom(c);
    Metadata.mesh = mesh(c);
    Metadata.paths = paths(c);
    
    save(paths(c).metafile,'Metadata');
end

return

stoopo
%%
figure(463)
plot([geom(succeeded).AR1],[geom(succeeded).AR2],'s','markerfacecolor','none','markeredgecolor','r','markersize',8);
hold on
plot([geom(~succeeded).AR1],[geom(~succeeded).AR2],'o','markerfacecolor','b','markeredgecolor','b','markersize',8);
hold off
grid on
%axis equal
axis tight

%%
figure(34)
inds = fliplr(1:590);
for i = 27:-1:1
    % if succeeded(i)
    meshname = paths(i).mesh_file;
    [Mesh, Metadata] = load_mesh(meshname);  %all paths should be the same for a given shape
    if isempty(Mesh)
        disp([meshname,'     didn''t work']);
        pause
        continue
    end
    figure(34)
    cla
    [s,e] = plot_mesh(Mesh);
    set(s,'facealpha',1);
    set(e,'edgealpha',0.3);
    title({paths(i).mesh_file,num2str(i),[num2str(Mesh.n_elem),'  elements']},'interpreter','none');
    grid off
    drawnow
    %     pause(0.25)
    
    % pause
    %didn't work:  4 12 15 16 17 21
    % 20 23
    %end
end