function [ A, A_motor_torque,A_force,A_torque,debug_info] = matrix_assembly_mex_wrapper( Mesh,matrix_props,assembly_input )
%wraps mexed function and assembles full A in native Matlab to avoid
%ridiculous limit for matrix size in Coder

%make sure we can even initialize the full damn thing before going further
A = zeros(matrix_props.n_rows, matrix_props.n_cols); %includes extra variables and equations for freeswim case


switch assembly_input.problemtype
    case 'forced'
        [Ax,Ay,Az,~,~,~,A_force,A_torque,debug_info] = matrix_assembly_mexed(Mesh,matrix_props,assembly_input);
        A(1:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Ax;
        %Ax = [];
        A(2:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Ay;
        %Ay = [];
        A(3:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Az;
        %Az = [];
        
        A = assembly_input.constants.multfactor*A;  %multiply by 1/8/pi as per eq. 3.4, DJ Smith "A boundary element regularized..." paper
        
        A_motor_torque = [];
    case 'freeswim'
        [Ax,Ay,Az,A_freeswim_rows,A_freeswim_cols,A_motor_torque,A_force,A_torque,debug_info] = matrix_assembly_mexed(Mesh,matrix_props,assembly_input);
        A(1:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Ax;
        %Ax = [];
        A(2:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Ay;
        %Ay = [];
        A(3:3:matrix_props.n_col*3,1:matrix_props.n_col*3) = Az;
        
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
            case 'dino'
                A(matrix_props.n_col*3 + 1:matrix_props.n_col*3 + 6, :) = A_freeswim_rows;
                A(:, matrix_props.n_col*3 + 1 : matrix_props.n_col*3 + 6) = A_freeswim_cols;
        end
        
        
end



