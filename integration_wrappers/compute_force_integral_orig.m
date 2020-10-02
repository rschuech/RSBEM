function [VL_unity, A_force] = compute_force_integral(Mesh, matrix_props, input)

%A_force is 3 x n_col*3 x length(Mesh)
%for total force on entire Mesh, use sum(A_force,3)

if length(Mesh) > 1 %multiple submeshes input, use global indices
    tot_elems = Mesh(end).indices.glob.unq_bounds.elem(2 );
else %only one submesh input, use local indices
    tot_elems = Mesh(1).n_elem;
end
%disp('right before integrate_unity');
[VL_unity] = integrate_unity(Mesh, input,'force');
%disp('right after integrate unity');

if nargout == 2
    %assemble matrix that can later be multiplied by the solution f to yield total
    %force on bug
    A_force = zeros(3,matrix_props.n_col*3,length(Mesh));
    %   disp('line 19');
    [~, elem_range] = bounds_from_global_inds(Mesh);
    
    [i_mesh_elems, ~] = global2local(1:tot_elems, Mesh, elem_range, 'elem');    %compute i_mesh_elem and local_elem, indices to current element's submesh and local element index
    
    %  disp('line 21');
    for i_mesh = 1:length(Mesh) %keep integrals for each submesh separate
        A_force_temp = zeros(3,matrix_props.n_col*3);
        %  disp('line 24');
        % parfor(elem_i = 1:tot_elems, input.performance.nthreads)  %hell, make it parallel anyway - hopefully worth it for huge runs
        %for some unknown reason, mex of the parfor causes a crash when run.
        % it gets all the way through the i_mesh loop but then crashes right
        % at the end
        for elem_i = 1:tot_elems  %not parallel but who cares - only need to do this once per run, and these are dumb operations
            %     disp('line 27');
            
            i_mesh_elem = i_mesh_elems(elem_i);
%             local_elem = local_elems(elem_i);  apparently unused
            
            
%             if length(Mesh) > 1 %multiple submeshes
% %                 [i_mesh_elem, local_elem] = global2local(elem_i, Mesh, elem_range,'elem');
% 
%             else
%                 i_mesh_elem = 1;
%                 local_elem = elem_i; %we're already looping over local elems in case of one submesh
%             end
            
            if i_mesh_elem ~= i_mesh %wrong submesh
                continue
            end
            
            col_inds = matrix_props.Col_inds(elem_i,1:18);
            % sum of forces = 0, entire object
            VL = VL_unity(elem_i,:);
            %disp('line 46');
            %VL is a 6 element vector that is the integral of hs * phi for this element
            %         A_force(1, col_inds(1:6))   = A_force(1, col_inds(1:6))   + VL;  %x, y, z components have same coeffs    x
            %         A_force(2, col_inds(7:12))  = A_force(2, col_inds(7:12))  + VL;  %x, y, z components have same coeffs    y
            %         A_force(3, col_inds(13:18)) = A_force(3, col_inds(13:18)) + VL;  %x, y, z components have same coeffs    z
            
            
            
            
            temp = zeros(3,matrix_props.n_col*3);
            temp(1,col_inds(1:6)) = VL; temp(2,col_inds(7:12)) = VL; temp(3,col_inds(13:18)) = VL;
            %disp('line 57');
            A_force_temp = A_force_temp + temp;
            %disp('line 59');
        end %elems
        
        A_force(:,:,i_mesh) = A_force_temp;
        %         disp('line 63');
    end
    %     disp('i_mesh loop end')
end
% disp('nargout if end')

