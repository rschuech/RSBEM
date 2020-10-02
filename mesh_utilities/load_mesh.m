function [Mesh, Metadata, original_node_indices] = load_mesh(meshname, varargin)
%loads mesh data that was created by Salome using gen_mesh.m
%specific to 6-node curved T6 triangles
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
    
    metaname = [meshname{1}(1:end-4),'_metadata.mat']; %metadata filename should always correspond to mesh filename
    
    if exist(metaname,'file')
        load(metaname);  %contains one variable, Metadata
        if isfield(Metadata,'mesh') && isfield(Metadata.mesh,'topology') && ~strcmp(Metadata.mesh.topology, 'sphere')
            disp(['Mesh topologically equivalent to ',Metadata.mesh.topology]);
        end
    else
%         if input.performance.verbose
        disp('Note:  Metadata file doesn''t exist');
%         end
%         Mesh = [];
        Metadata = [];
%         return
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
    input.performance.randomize_nodes = true;
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

Mesh.n_nodes = temp(1,1);  %# nodes, may or may not be the number that's actually used - see below
Mesh.n_elements = temp(1,2); %# elements

original_node_indices = temp(2:Mesh.n_nodes+1,1);
% original_indices.elem = temp(Mesh.n_nodes+2:end,1);

Mesh.nodes = temp(2:Mesh.n_nodes+1,2:4); %node coords

Mesh.elements = temp(Mesh.n_nodes+2:end,3:end);  %node indices that form each element

%remove non-curved-triangular elements with less than 6 nodes
inds = sum(Mesh.elements == 0,2) > 0;
% original_indices.elem(inds) = [];
Mesh.elements(inds,:) = [];  %some elements may be degenerate, not having 6 nodes - remove them
Mesh.n_elements = size(Mesh.elements,1); %this is the actual # elements

node_inds_used = unique(Mesh.elements(:));  %sometimes not all listed nodes are
%actually used in any elements, so these are the ones actually used...
[~,inds] = ismember(node_inds_used, original_node_indices);  %row inds of original node_inds list corresponding to node_inds actually used
Mesh.nodes = Mesh.nodes(inds(inds ~= 0),:);  %keep only the nodes actually used in elements
original_node_indices = original_node_indices(inds(inds ~= 0));
Mesh.n_nodes = size(Mesh.nodes,1); %actual # nodes

%now fix node indices defining each element to correspond to new row indices
%created when Mesh.nodes was modified above

% newnum = 0; %index matching row of new Mesh.nodes
% for nodenum = nodenums' %step through old, unique node inds
%     newnum = newnum + 1;  %new row index of current node
%     Mesh.elements(Mesh.elements == nodenum) = newnum; %update node ind wherever it's mentioned in an element
% end
%don't need above anymore since we are keeping node_inds and elem_inds
%around instead of using the matrix row numbers

%load_mesh_wrapper(meshname, input, persistant_data(si).rand_inds, input_diffs); 
if input.performance.randomize_nodes
    if nargin >= 3 && ~isempty(varargin{2})
        rand_inds = varargin{2};
       
    else
        rng(1);  % always generate same randomization of a given mesh for reproducibility
        rand_inds = randperm(Mesh.n_nodes);
    end
    [~, inds] = sort(rand_inds);  %map from orig to new randomized row inds
    Mesh.nodes = Mesh.nodes(rand_inds,:); %randomize node order
    original_node_indices = original_node_indices(rand_inds);
    % Mesh.elements = inds(Mesh.elements);  %update elements to use new node inds
    %     don't need to do anymore since we have Mesh.node_inds(?)
    Metadata.mesh.rand_inds = rand_inds;
end


Mesh.nodes = Mesh.nodes * rescale;  % rescale units of geometry if desired

% In the future we may split up closed surfaces into multiple submeshes and there is no easy way for the code to determine that the individual
% submeshes make up a closed surface (e.g., what if dino body were split in two, but transverse/hairs were still attached).  Since mesh topology isn't
% really used except in deciding whether to incl DL integrals, just don't bother talking about submesh topology and instead rely on hardcoded eliminate_DL
% flags in settings_inputs, trusting user to use knowledge of open/closed topologies in setting those.


% edges = get_edges(Mesh);
% n_neighbors = cellfun(@length,edges.elements);  % number of elements each edge is a member of
% types = unique(n_neighbors);
% if isequal(types,2) % each edge belongs to 2 elements
%     Mesh.topology = 'closed';
% elseif isqual(types,[1 2]) % some edges belong to 1 element, must be the edge of an open surface
%     Mesh.topology = 'open';
% else
%     error('Mesh topology is weird');
% end

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

