function out = integrand(X,integrand_constants)

%#codegen

% all this could be sent in as a function handle to adsimp instead of
% buried here?

% varargin is for use with integral_type = 'reg_stresslet':
% varargin{1} is type of reg_stresslet integrands, i.e. {'Tn' , 'rTn' , 'phiTn} or a subset (but with this order preserved in output)
% varargin{2} is whether to use singularity-removal identity "trick" for the reg_stresslet integrals (true or false)
% this identity allows one to integrate (u - u_c)Tn [instead of uTn] where u_c is velocity at collocation pt

xi = X(1);  eta = X(2);
out = [];

switch integrand_constants.integral_type
    
    case "reg_stresslet"
        
        [x_interp, phi, hS, n] = T6interp(integrand_constants.element_nodes,xi,eta,integrand_constants.shape_parameters);
        T = calcT(x_interp,integrand_constants.x_col,integrand_constants.eps2);
%         r = x_interp - integrand_constants.refpoint;
        
        %         if any( ismember({'Tn','rTn','phiTn'},varargin{1}) )
        
%         if isequal( varargin{1}{:}, 'Tn') && varargin{2} % trick and we're trying to integrate UTn *only* where U is a constant translational velocity
%             Tn = zeros(3,3); % the integral of Tn is not actually zero, but the integral of (U - U_c)Tn will be zero since U is constant throughout the
            % geometry, so we can set Tn = 0 here to achieve this end result
%         else % only want to integrate UTn with no trick, or, we first need Tn to integrate rTn or phiTn
            Tn = NaN(3,3);
            %             for i = 1:3
            %                 for j = 1:3
            %                     Tn(i,j) = sum(squeeze(T(i,j,:)) .* n);  % 3 x 3
            %                 end
            %             end
            
            % need to check if this is faster (mexed) than above?
            for l = [1 2 3 5 6 9] % linear indices of lower triangular entries
                [i,j] = ind2sub([3,3],l);
                Tn(l) = sum(squeeze(T(i,j,:)) .* n);  % 3 x 3
            end
            Tn(4) = Tn(2); Tn(7) = Tn(3);  Tn(8) = Tn(6); % Tn is symmetric
            
            % or, would writing out each of 9 entries be fastest: ??
            %             Tn(1,1) = sum(squeeze(T(1,1,:)) .* n);
            %             Tn(2,1) = sum(squeeze(T(2,1,:)) .* n);
            %             Tn(3,1) = sum(squeeze(T(3,1,:)) .* n);
            %             Tn(1,2) = Tn(2,1);
            %             Tn(2,2) = sum(squeeze(T(2,2,:)) .* n);
            %             Tn(3,2) = sum(squeeze(T(3,2,:)) .* n);
            %             Tn(1,3) = Tn(3,1);
            %             Tn(2,3) = Tn(3,2);
            %             Tn(3,3) = sum(squeeze(T(3,3,:)) .* n);
            %
%         end % no Tn = [] option because if we don't want Tn, we must want either rTn or phiTn or both, which require Tn as an intermediate calc
        % Tn are the coeffs that multiply U when computing int(U*T*n) or that multiply unknown u at each collocation pt for free-slip problem
        % only need Tn if we are not using singularity removal identity, since this term (U - U_c)Tn cancels out otherwise since U = U_c everywhere
        %         end
        
        
        % We don't ever actually need uTn since whether u is known or unknown, we can factor it out by writing u as sum(phi.*u_v) or with trick,
        % sum(phi.*(u_v - u_c)), and integrate phiTn instead.  For the trick, we can hardcode out the hardest integrals since we know at those verts, the
        % u_v = u_c and the integral for the corresponding phi term will be zero
        % factoring u out is advantageous since the scale of u is taken out of the integrals and we can base adaptive integration error tolerances solely
        % on the scale of the geometry
        
        % uTn are the constants that go on the RHS when computing
        % int(u*T*n) where u is a known velocity at x_interp, either
        % relative flagella velocity in mobility problem or imposed
        % body velocity in resistance problem
        
        %         if ismember('uTn',varargin{1})
        %             u_interp = integrand_constants.velocities' * phi';  % (3 x 6)  *  (6 x 1) = (3 x 1)
        %
        %             if varargin{2} % use singularity removal so we integrate (u - u_c)Tn
        %                 u_mod = u_interp - integrand_constants.u_col;
        %             else
        %                 u_mod = u_interp;
        %             end
        %
        %             uTn = (u_mod' * Tn)';  % 3 x 1
        %         else
        %             uTn = [];
        %         end
        
        
%         if ismember('rTn',varargin{1})
%             % rTn are the coeffs that multiply Omega when computing
%             %          int((Omega x r)*T*n)
%             
%             if varargin{2} % use singularity removal so that (Omega x r - Omega x r_c)Tn = ( Omega x (r - r_c) )Tn
%                 r_c = integrand_constants.x_col - integrand_constants.refpoint;
%                 r_mod = r - r_c;
%             else
%                 r_mod = r;
%             end
%             rTn = [r_mod(2)*Tn(1,3) - r_mod(3)*Tn(1,2) , r_mod(3)*Tn(1,1) - r_mod(1)*Tn(1,3) , r_mod(1)*Tn(1,2) - r_mod(2)*Tn(1,1);...
%                 r_mod(2)*Tn(2,3) - r_mod(3)*Tn(2,2) , r_mod(3)*Tn(2,1) - r_mod(1)*Tn(2,3) , r_mod(1)*Tn(2,2) - r_mod(2)*Tn(2,1);...
%                 r_mod(2)*Tn(3,3) - r_mod(3)*Tn(3,2) , r_mod(3)*Tn(3,1) - r_mod(1)*Tn(3,3) , r_mod(1)*Tn(3,2) - r_mod(2)*Tn(3,1); ];  % 3 x 3
%         else
%             rTn = [];
%         end
        
%         if ismember('phiTn',varargin{1}) % for free slip BC when u is unknown and must be solved for
            phi_mod = phi;
            if integrand_constants.DL_singularity_removal % with trick, if this elem contains the collocation pt, the integrals involving (u_v - u_c) at that vertex will be zero,
                % regardless of phi for that vertex.  A shortcut to zero these integrals is to set phi to zero for them even though the zero doesn't come
                % from the phi.  col_vert_ind will either be the local index (1 - 6) of the collocation vert, or [] if this element doesn't contain the
                % collocation pt, in which case phi_mod = phi.
                phi_mod(integrand_constants.coll_node_ind) = 0;
            end
            
            phiTn = kron(Tn(:) , phi_mod); % 54 x 1 due to 6 phi and 9 Tn entries   kron(9 x 1 , 6 x 1) = 54 x 1
            % phi.*Tn(1,1); phi.*Tn(2,1); phi.*Tn(3,1); phi.*Tn(1,2); phi.*Tn(2,2) .... phi.*Tn(3,3)
            % think of moving across rows of uTn with u = sum(phi.*u_v) substituted, u_v are velocities at elem verts
            % again we are repeating indentical integrals here due to symmetry of Tn (6 unique entries) - there are 6*6 = 36 unique integrals out of 54 total
            % but not repeating integrals didn't seem to improve speed with the SL reg Stokeslet below, so perhaps wouldn't help here either?
%         else
%             phiTn = [];
%         end
        
        
        
%         if ~ismember('Tn',varargin{1}) % even if we calculated Tn in order to calculate rTn or phiTn, delete it here if we don't want to integrate Tn
%             Tn = [];
%         end
        
%         out = hS*[Tn(:); rTn(:); phiTn];
        out = hS*phiTn;
        
        
        
    case "reg_stokeslet"
        [x_interp, phi, hS] = T6interp(integrand_constants.element_nodes,xi,eta,integrand_constants.shape_parameters);
        S = calcS(x_interp,integrand_constants.x_col,integrand_constants.eps2);
        
        out =  hS*[S(1,1)*phi; S(1,2)*phi; S(1,3)*phi; ...
            S(2,1)*phi; S(2,2)*phi; S(2,3)*phi; ...
            S(3,1)*phi; S(3,2)*phi; S(3,3)*phi];
        
        %above wastes time since S is symmetric and S(2,1), S(3,1), S(3,2) are repeated, so the
        %integrals are repeated    ....  and yet, the mexed version is still
        %faster than a mexed version that copies these values from S(1,2),
        %S(1,3), S(2,3)
        
        % check this again?
        
        %      out =  hS*  [S(1,1)*phi S(1,2)*phi S(1,3)*phi ...
        %                              S(2,2)*phi S(2,3)*phi ...
        %                                         S(3,3)*phi]';
%     case "S_dot_f"  % never used?
%         % integrates S dot f or equivalently f dot S (S is symmetric so
%         % only difference is whether f and the result is column or row
%         % vector) to compute velocity at any field points
%         
%         % it appears this is currently used for field velocity calc?
%         % but it would be more logical to use the reg_Stokeslet integrand for field velocity also - should reuse matrix_assembly_mex but instead of
%         % building a matrix, can compute one (or 3?) rows at a time for each field point, and immediately multiply by known f to obtain the field velocity
%         % at each field point.  No point in saving a huge matrix here since it probably can't be reused at all?  (Since the field point is never on the
%         % mesh, we can't ever use rigid rotation trick to skip any integrals?)
%         
%         [x_interp, phi, hS] = T6interp(integrand_constants.element_nodes,xi,eta,integrand_constants.shape_parameters);
%         % [f_interp] = T6interp(integrand_constants.f,xi,eta,shape_parameters);
%         f_interp = integrand_constants.f' * phi;
%         
%         S = calcS(x_interp,integrand_constants.x_field,integrand_constants.eps2);
%         
%         out =  hS*[S(1,1) S(1,2) S(1,3) ; ...
%             S(2,1) S(2,2) S(2,3) ; ...
%             S(3,1) S(3,2) S(3,3)] * f_interp;
        
        
    case "unity"
        [~,  phi, hS] = T6interp(integrand_constants.element_nodes,xi,eta,integrand_constants.shape_parameters);
        out =  hS*phi;  %only 6 numbers - same coeffs are repeated for integrating fx, fy, and fz force components in the 3 component equations for the force balance
    case "moments"
        [x_interp, phi, hS] = T6interp(integrand_constants.element_nodes,xi,eta,integrand_constants.shape_parameters);
        r = x_interp - integrand_constants.refpoint;  %want r vector measured from point along center axis, refpoint
        
        out =  hS*[r(1) * phi;    r(2) * phi;    r(3) * phi];
        
        % cross product x_interp cross F involves 3 components of x_interp
        % multiplied in various permutations with the 3 components of F, each
        % coming from phi.  So only need to integrate the above values, and can
        % later insert these coeffs at the appropriate places in A to go with the
        % correct force components
        
    case "flux" %use divergence theorem to compute mesh volume from a surface integral
        % integral of F dot n * hS * phi
        %for now, only used to calculate volume, so hardcoded F = [1/3 x   1/3 y  1/3 z]
        [x_interp, phi, hS, n] = T6interp(integrand_constants.element_nodes,xi,eta,integrand_constants.shape_parameters);
        out = 1/3*(x_interp' * n) * hS * phi;  % (1 x 3) * (3 x 1) * (1 x 1) * (6 x 1) = (6 x 1)
        
    case "centroid"  %uses divergence theorem with F = [1/2 x^2  0   0] or F = [0  1/2 y^2  0] or F = [0  0  1/2 z^2]
        %to calculate integrals for centroid of volume enclosed by mesh
        %assuming uniform unit density
        [x_interp, phi, hS, n] = T6interp(integrand_constants.element_nodes,xi,eta,integrand_constants.shape_parameters);
        out = [  [1/2*x_interp(1)^2  0  0] * n * hS * phi; ...
            [0  1/2*x_interp(2)^2  0] * n * hS * phi; ...
            [0  0  1/2*x_interp(3)^2] * n * hS * phi;   ];
        
    case "signed_flow" % similar to flux but need to use computed velocity instead of simple hardcoded F
        % because we care about sign of V dot n here, we can't output a
        % matrix and reuse it for different V, but need to compute V dot n
        % right here (at expense of needing to re-integrate for different
        % phases)
        [x_interp, phi, hS, n] = T6interp(integrand_constants.element_nodes,xi,eta,integrand_constants.shape_parameters);
        % assuming field_vel is 6 x 3
        
        V = sum( integrand_constants.field_vel .* repmat(phi,1,3) , 1 );   % sum --> 1 x 3   V at current xi, eta point
        % n should be 3 x 1
        local_flow = V * n * hS;
        
        %        out = [local_flow local_flow]';
        
        if local_flow > 0
            out =  [local_flow 0]';
        else
            out = [0 local_flow]';
        end
        
end