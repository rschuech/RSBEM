
function [Mesh] = store_mesh_constants_mex(Mesh,input)

for i = 1:length(Mesh)  %loop over submeshes in case there is more than one
    
    %all integrals below are gratuitously adaptive, using tolerances in input struct
    
    centroid = NaN(3,size(Mesh(i).elems,1)); %centroid of each element
    
    %     Mesh(i).area = NaN(0,1);  %area of each element
    elem_params = NaN(3,size(Mesh(i).elems,1));
    %     Mesh(i).Centroid = NaN(3,1);  %centroid of entire mesh volume assuming uniform density
    
    parfor (elem_i = 1:size(Mesh(i).elems,1), input.performance.nthreads)
        %for elem_i = 1:size(Mesh.elems,1)
        %elem_i / size(elems,1)
        subverts = Mesh(i).verts(Mesh(i).elems(elem_i,:),:);
        
        alpha = 1/(1+norm(subverts(4,:)-subverts(2,:))/norm(subverts(4,:)-subverts(1,:)));
        beta = 1/(1+norm(subverts(6,:)-subverts(3,:))/norm(subverts(6,:)-subverts(1,:)));
        gamma = 1/(1+norm(subverts(5,:)-subverts(2,:))/norm(subverts(5,:)-subverts(3,:)));
        
        elem_params(:,elem_i) = [alpha beta gamma];
        
        centroid(:,elem_i) = T6interp(subverts,1/3,1/3,elem_params(:,elem_i));
        
        
    end
    Mesh(i).elem_params = elem_params;
    Mesh(i).centroid = centroid;
    
  %  i
    Mesh(i).area = ones(Mesh(i).n_elem,1); %need to temporarily initialize this as something since abstol for integration is always scaled by element area, but we're trying to calculate element area....
   %'starting area'
    [VL_unity] = integrate_unity(Mesh(i), input,'area'); %compute integrals of hS * phi
%    'area done'
    Mesh(i).area = sum(VL_unity,2);  %area is just integral of 1 * hS * phi so can just add up all the hS * phi integrals for each element
%     'starting flux'
    VL_flux = integrate_flux(Mesh(i), input, 'volume');  %compute integrals of F dot n * hS * phi ala divergence theorem, F hardcoded as [1/3 x  1/3 y  1/3 z] but other F can be added later
%     'flux done'
    Mesh(i).Volume = sum(VL_flux(:));  %not sure integrals for each individual element are very useful here, but the whole thing is the total volume
    
%     'starting centroid'
    [VL_centroid] = centroid_integrals(Mesh(i), input);  %output:  row is elem, cols are phi = 1..6 for x, y, z centroid integrals
%    'centroid done'
    Mesh(i).Centroid(1,1) = sum(sum(VL_centroid(:,1:6))) / Mesh(i).Volume; %not sure if finer grained sums mean anything, but complete sum over all elements is integral (x dV)  and then normalize by total volume
    Mesh(i).Centroid(2,1) = sum(sum(VL_centroid(:,7:12))) / Mesh(i).Volume;
    Mesh(i).Centroid(3,1) = sum(sum(VL_centroid(:,13:18))) / Mesh(i).Volume;
    
    
    % Mesh_all_out(i) = Mesh; %put back into combined Mesh_all variable for output
    
end