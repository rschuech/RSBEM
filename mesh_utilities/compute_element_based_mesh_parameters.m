
function [Mesh] = compute_element_based_mesh_parameters(Mesh, input)

% computes parameters corresponding to each element (e.g. centroid, area) as well as a few properties of the entire mesh (e.g. centroid, volume) and stores
% them in Mesh


%#codegen

for i = 1:length(Mesh)  %loop over submeshes in case there is more than one
    
    %all integrals below are gratuitously adaptive, using tolerances in input struct
    
    centroid = NaN(3,size(Mesh(i).elements,1)); %centroid of each element
    
    shape_parameters = NaN(size(Mesh(i).elements,1),3);
    
    
    %     parfor (elem_i = 1:size(Mesh(i).elements,1), input.performance.nthreads)
    for elem_i = 1:size(Mesh(i).elements,1) %parfor broke when normals and tangents calc was added, prolly not worth fixing
        %elem_i / size(elements,1)
        node_coords = Mesh(i).nodes(Mesh(i).elements(elem_i,:),:); % coords of the 6 nodes in this element
        % see Pozrikidis page 123 for what alpha, beta, gamma mean:
        alpha = 1/(1+norm(node_coords(4,:)-node_coords(2,:))/norm(node_coords(4,:)-node_coords(1,:))); %corresponds to location of edge node # 4
        beta = 1/(1+norm(node_coords(6,:)-node_coords(3,:))/norm(node_coords(6,:)-node_coords(1,:))); % corresponds to location of edge node # 6
        gamma = 1/(1+norm(node_coords(5,:)-node_coords(2,:))/norm(node_coords(5,:)-node_coords(3,:))); % corresponds to location of edge node # 5
        
        shape_parameters(elem_i,:) = [alpha beta gamma];
        
        centroid(:,elem_i) = T6interp(node_coords,1/3,1/3,shape_parameters(elem_i,:)); %coords of element centroid
        
    end % elements
    
    Mesh(i).shape_parameters = shape_parameters;
    Mesh(i).centroid = centroid;
    
    Mesh(i).area = ones(Mesh(i).n_elements,1); %need to temporarily initialize this as something since abstol for integration is always scaled by element area, but we're trying to calculate element area....
    %'starting area'
    %     [VL_unity] = integrate_unity(Mesh(i), input,'area'); %compute integrals of hS * phi
    %     'starting flux'
    %     VL_flux = integrate_flux(Mesh(i), input, 'volume');  %compute integrals of F dot n * hS * phi ala divergence theorem, F hardcoded as [1/3 x  1/3 y  1/3 z] but other F can be added later
    %     'flux done'
    %     Mesh(i).Volume = sum(VL_flux(:));  %not sure integrals for each individual element are very useful here, but the whole thing is the total volume
    
    %     'starting centroid'
    %     [VL_centroid] = centroid_integrals(Mesh(i), input);  %output:  row is elem, cols are phi = 1..6 for x, y, z centroid integrals
    %     %    'centroid done'
    %     Mesh(i).Centroid(1,1) = sum(sum(VL_centroid(:,1:6))) / Mesh(i).Volume; %not sure if finer grained sums mean anything, but complete sum over all elements is integral (x dV)  and then normalize by total volume
    %     Mesh(i).Centroid(2,1) = sum(sum(VL_centroid(:,7:12))) / Mesh(i).Volume;
    %     Mesh(i).Centroid(3,1) = sum(sum(VL_centroid(:,13:18))) / Mesh(i).Volume;
    
end % submeshes

[VL_unity] = integration_wrapper(Mesh, 'unity', ...
    input.accuracy.integration_tol.area.abstol,...
    input.accuracy.integration_tol.area.reltol,...
    input.accuracy.integration_tol.area.maxevals,...
    input.performance.nthreads); %compute integrals of hS * phi
[VL_flux] = integration_wrapper(Mesh, 'flux', ...
    input.accuracy.integration_tol.volume.abstol,...
    input.accuracy.integration_tol.volume.reltol,...
    input.accuracy.integration_tol.volume.maxevals,...
    input.performance.nthreads); %compute integrals of F dot n * hS * phi ala divergence theorem, F hardcoded as [1/3 x  1/3 y  1/3 z] but other F can be added later
[VL_centroid] = integration_wrapper(Mesh, 'centroid', ...
    input.accuracy.integration_tol.centroid.abstol,...
    input.accuracy.integration_tol.centroid.reltol,...
    input.accuracy.integration_tol.centroid.maxevals,...
    input.performance.nthreads); %output:  row is elem, cols are phi = 1..6 for x, y, z centroid integrals

for i = 1:length(Mesh)
    Mesh(i).area = sum(VL_unity{i},2);  %area is just integral of 1 * hS * phi so can just add up all the hS * phi integrals for each element
    Mesh(i).Volume = sum(VL_flux{i}(:)); %not sure integrals for each individual element would ever be useful here, but the whole thing is the total volume
    for j = 1:3
        Mesh(i).Centroid(j,1) = sum(sum(VL_centroid{i}(:,6*(j-1)+1:6*(j-1)+6))) / Mesh(i).Volume; %not sure if finer grained sums mean anything, but complete sum over all elements is integral (x dV)  and then normalize by total volume
    end
end