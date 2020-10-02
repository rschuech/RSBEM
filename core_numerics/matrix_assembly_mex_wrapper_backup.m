function [ A, A_motor_torque,A_force,A_torque,debug_info] = matrix_assembly_mex_wrapper_backup( Mesh,matrix_props,assembly_input )
%wraps mexed function and assembles full A in native Matlab to avoid
%ridiculous limit for matrix size in Coder

%make sure we can even initialize the full damn thing before going further
A = zeros(matrix_props.n_rows, matrix_props.n_cols); %includes extra variables and equations for freeswim case


switch assembly_input.problemtype
    case 'forced'
        [Ax_1,Ay_1,Az_1,  Ax_2,Ay_2,Az_2,  Ax_3,Ay_3,Az_3,  Ax_4,Ay_4,Az_4,       ~,~,~,A_force,A_torque,debug_info] = matrix_assembly_mexed(Mesh,matrix_props,assembly_input);
        A(1:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = [Ax_1 Ax_2 Ax_3 Ax_4];  %Ax is split up further over columns
        %Ax = [];
        A(2:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = [Ay_1 Ay_2 Ay_3 Ay_4];
        %Ay = [];
        A(3:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = [Az_1 Az_2 Az_3 Az_4];
        %Az = [];
        
        A = assembly_input.constants.multfactor*A;  %multiply by 1/8/pi as per eq. 3.4, DJ Smith "A boundary element regularized..." paper
        
        A_motor_torque = [];
    case 'freeswim'
        [Ax_1,Ay_1,Az_1,  Ax_2,Ay_2,Az_2,  Ax_3,Ay_3,Az_3,  Ax_4,Ay_4,Az_4,       A_freeswim_rows,A_freeswim_cols,A_motor_torque,A_force,A_torque,debug_info] = matrix_assembly_mexed(Mesh,matrix_props,assembly_input);
        A(1:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = [Ax_1 Ax_2 Ax_3 Ax_4];
        %Ax = [];
        A(2:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = [Ay_1 Ay_2 Ay_3 Ay_4];
        %Ay = [];
        A(3:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = [Az_1 Az_2 Az_3 Az_4];
        
        A = assembly_input.constants.multfactor*A;  %multiply by 1/8/pi as per eq. 3.4, DJ Smith "A boundary element regularized..." paper
        
        switch assembly_input.bugtype
            case 'bacteria'
                switch assembly_input.tail.motorBC
                    case 'freq'
                        A(matrix_props.n_col*3 + 1:matrix_props.n_col*3 + 6, :) = A_freeswim_rows;
                        A(:, matrix_props.n_col*3 + 1 : matrix_props.n_col*3 + 6) = A_freeswim_cols;
                    case 'torque'
                        A(matrix_props.n_col*3 + 1:matrix_props.n_col*3 + 7, :) = A_freeswim_rows;
                        A(:, matrix_props.n_col*3 + 1 : matrix_props.n_col*3 + 7) = A_freeswim_cols;
                end
            case {'dino', 'sheet'}
                A(matrix_props.n_col*3 + 1:matrix_props.n_col*3 + 6, :) = A_freeswim_rows;
                A(:, matrix_props.n_col*3 + 1 : matrix_props.n_col*3 + 6) = A_freeswim_cols;
        end
        
        
end



