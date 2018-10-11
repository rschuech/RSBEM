function [Mesh, Metadata] = load_mesh(meshname, varargin)
%loads mesh data that was created by Salome using gen_mesh.m
%specific to 6-vertex curved T6 triangles
%also loads metadata file with geometry and mesh refinement info

%takes input struct as an optional 2nd input; defaults are used if omitted

%rand_inds is optional 3rd input - if a forced or freeswim run was just
%done and we're doing the other now, use existing rand_inds so that most of
%A matrix can be copied instantly

%4th optional input is either 'metadata', which if true, causes us to only
%attempt to load Metadata and then return, or 'mesh', which causes us to
%only attemt to load Mesh and then return


%IMPORTANT:  It's best not to try to work with individual submeshes since
%many things, like the integrals in store_mesh_constants, might break.
%Instead, if interested in one or more submeshes, reload just these into a
%new Mesh variable and work with that.

if ~( nargin >= 4 && strcmp(varargin{3}, 'mesh'))  %we've NOT input Mesh-only flag, don't skip metadata
    
    metaname = [meshname(1:end-4),'_metadata.mat']; %metadata filename should always correspond to mesh filename
    
    if exist(metaname,'file')
        load(metaname);  %contains one variable, Metadata
        if isfield(Metadata.mesh,'topology') && ~strcmp(Metadata.mesh.topology, 'sphere')
            disp(['Mesh topologically equivalent to ',Metadata.mesh.topology]);
        end
    else
%         if input.performance.verbose
%         disp('Metafile doesn''t exist; aborting load_mesh');
%         end
        Mesh = [];
        Metadata = [];
        return
    end
    
    
    if nargin >= 4 && strcmp(varargin{3}, 'metadata')  %we've input Metadata_only flag
        Mesh = [];
        return
    end
    
    
else %skip metadata
    
    Metadata = [];  %skipped, only try to load Mesh
    
end


if ~exist(meshname,'file')
%     if input.performance.verbose
%     disp('Mesh file doesn''t exist');
%     end
    Mesh = [];
    return
end

if nargin >= 2 && ~isempty(varargin{1})
    input = varargin{1};
else %use some defaults
    input.performance.randomize_verts = true;
    input.performance.nthreads = feature('numCores');
    %     input.accuracy.integration_tol.area.reltol = 1E-1;
    %     input.accuracy.integration_tol.volume.reltol = 1E-1;
    %     input.accuracy.integration_tol.centroid.reltol = 1E-1;
    names = {'area','volume','centroid'};
    for n = 1:length(names)
        input.accuracy.integration_tol.(names{n}).abstol = 0;
        input.accuracy.integration_tol.(names{n}).reltol = 1E-1;
        input.accuracy.integration_tol.(names{n}).maxevals = Inf;
    end
end

try
rescale = 1; %rescale geometry by this factor
%current convention is to use microns for all lengths
temp = dlmread(meshname);

Mesh.n_vert = temp(1,1);  %# vertices, may or may not be the number that's actually used - see below
Mesh.n_elem = temp(1,2); %# elements

Mesh.indices.orig.vert = temp(2:Mesh.n_vert+1,1);
Mesh.indices.orig.elem = temp(Mesh.n_vert+2:end,1);

Mesh.verts = temp(2:Mesh.n_vert+1,2:4); %vertex coords

Mesh.elems = temp(Mesh.n_vert+2:end,3:end);  %vertex indices that form each element

%remove non-curved-triangular elements with less than 6 nodes
inds = sum(Mesh.elems == 0,2) > 0;
Mesh.indices.orig.elem(inds) = [];
Mesh.elems(inds,:) = [];  %some elements may be degenerate, not having 6 vertices - remove them
Mesh.n_elem = size(Mesh.elems,1); %this is the actual # elements

vertinds_used = unique(Mesh.elems(:));  %sometimes not all listed vertices are
%actually used in any elements, so these are the ones actually used...
[~,inds] = ismember(vertinds_used, Mesh.indices.orig.vert);  %row inds of original vert_inds list corresponding to vertinds actually used
Mesh.verts = Mesh.verts(inds(inds ~= 0),:);  %keep only the verts actually used in elems
Mesh.indices.orig.vert = Mesh.indices.orig.vert(inds(inds ~= 0));
Mesh.n_vert = size(Mesh.verts,1); %actual # vertices

%now fix vert indices defining each element to correspond to new row indices
%created when Mesh.verts was modified above

% newnum = 0; %index matching row of new Mesh.verts
% for vertnum = vertnums' %step through old, unique vert inds
%     newnum = newnum + 1;  %new row index of current vertex
%     Mesh.elems(Mesh.elems == vertnum) = newnum; %update vert ind wherever it's mentioned in an element
% end
%don't need above anymore since we are keeping vert_inds and elem_inds
%around instead of using the matrix row numbers

%load_mesh_wrapper(meshname, input, persistant_data(si).rand_inds, input_diffs); 
if input.performance.randomize_verts
    if nargin >= 3 && ~isempty(varargin{2})
        rand_inds = varargin{2};
       
    else
        rand_inds = randperm(Mesh.n_vert);
    end
    [~, inds] = sort(rand_inds);  %map from orig to new randomized row inds
    Mesh.verts = Mesh.verts(rand_inds,:); %randomize vert order
    Mesh.indices.orig.vert = Mesh.indices.orig.vert(rand_inds);
    % Mesh.elems = inds(Mesh.elems);  %update elements to use new vert inds
    %     don't need to do anymore since we have Mesh.vert_inds(?)
    Metadata.mesh.rand_inds = rand_inds;
end


Mesh.verts = Mesh.verts * rescale;  % rescale units of geometry if desired

%don't need anymore?  
%[Mesh] = global_inds(Mesh);  %here there is only one submesh that we just loaded, but create placeholders that will be updated later

% create smaller input struct with just necessary parameters to reduce mex
% pain
% temp_input.performance.nthreads = input.performance.nthreads;
% temp_input.accuracy.integration_tol.area = input.accuracy.integration_tol.area;
% temp_input.accuracy.integration_tol.centroid = input.accuracy.integration_tol.centroid;
% temp_input.accuracy.integration_tol.volume = input.accuracy.integration_tol.volume;

%[Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
catch
    disp('Mesh loading failed for unexpected reason')
    Mesh = [];
end

