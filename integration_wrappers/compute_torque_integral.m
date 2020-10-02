function [VL_moments, A_torque] = compute_torque_integral(Mesh, matrix_props, input)

%unlike other integration functions, this and compute_force_integral are
%designed to only take entire Mesh as input, since it must go with solution
%f to be useful and f is presumably for entire Mesh

%A_torque is 3 x n_col*3 x length(Mesh)
%for total torque on entire Mesh, use sum(A_torque,3)

if length(Mesh) > 1 %multiple submeshes input, use global indices
    tot_elems = Mesh(end).indices.glob.unq_bounds.elem(2 );
else %only one submesh input, use local indices
    tot_elems = Mesh(1).n_elem;
end

[VL_moments] = integrate_moments(Mesh, input,'torque');

if nargout == 2
    %assemble matrix that can later be multiplied by the solution f to yield total
    %force and total torque on bug
    A_torque = zeros(3,matrix_props.n_collocation*3,length(Mesh));
    [~, elem_range] = bounds_from_global_inds(Mesh);
     [i_mesh_elems, ~] = global2local(1:tot_elems, Mesh, elem_range, 'elem');    %compute i_mesh_elem and local_elem, indices to current element's submesh and local element index
    
    for i_mesh = 1:length(Mesh) %keep integrals for each submesh separate
        
        A_torque_temp = zeros(3,matrix_props.n_collocation*3);
        %parfor (elem_i = 1:tot_elems, input.performance.nthreads)  %hell, make it parallel anyway - hopefully worth it for huge runs
        for elem_i = 1:tot_elems  %not parallel but who cares - only need to do this once per run, and these are dumb operations
            
            i_mesh_elem = i_mesh_elems(elem_i);
          %  local_elem = local_elems(elem_i);  apparently unused
            
            %             if length(Mesh) > 1 %multiple submeshes
            %                 [i_mesh_elem, local_elem] = global2local(elem_i, Mesh, elem_range,'elem');
            %             else
            %                 i_mesh_elem = 1;
            %                 local_elem = elem_i; %we're already looping over local elems in case of one submesh
%             end
            
            if i_mesh_elem ~= i_mesh %wrong submesh
                continue
            end
            
            col_inds = matrix_props.Col_inds(elem_i,1:18);
            
            %sum of torques = 0, entire object (measured around refpoint, which moves with
            %object and is along motor axis)
            
            VL = VL_moments(elem_i,:);
            %VL is a 18 element vector that is the integral of hs * (x*phi  y*phi  z * phi) for this element
            
            %each component of r cross F contains two components of W and r, hence 2
            %lines per component of r cross F
            %         A_torque(1, col_inds(13:18)) = A_torque(1, col_inds(13:18)) + VL(7:12);
            %         A_torque(1, col_inds(7:12))  = A_torque(1, col_inds(7:12))  - VL(13:18);  %x component of xinterp cross F
            %         A_torque(2, col_inds(1:6))   = A_torque(2, col_inds(1:6))   + VL(13:18);
            %         A_torque(2, col_inds(13:18)) = A_torque(2, col_inds(13:18)) - VL(1:6);   %y component of xinterp cross F
            %         A_torque(3, col_inds(7:12))  = A_torque(3, col_inds(7:12))  + VL(1:6);
            %         A_torque(3, col_inds(1:6))   = A_torque(3, col_inds(1:6))   - VL(7:12);    %z component of xinterp cross F
            %
            
            temp = zeros(3, matrix_props.n_collocation*3);
            temp(1,col_inds(13:18)) = VL(7:12);
            temp(1,col_inds(7:12)) = -VL(13:18);
            temp(2,col_inds(1:6)) = VL(13:18);
            temp(2,col_inds(13:18)) = -VL(1:6);
            temp(3,col_inds(7:12)) = VL(1:6);
            temp(3,col_inds(1:6)) = -VL(7:12);
            
            A_torque_temp = A_torque_temp + temp;
            
            
        end %elems
        A_torque(:,:,i_mesh) = A_torque_temp;
        
    end
    
end