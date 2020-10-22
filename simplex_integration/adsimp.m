function [INTEGRALS, ABS_ERROR_ESTIMATE, NUM_F_EVALS_USED, FL] = adsimp( ndims, VRTS, ncomponents,  max_F_evals, ERROR_ABS, ERROR_REL, rule_constants, integrand_constants)

%#codegen
% coder.extrinsic('tic');
% coder.extrinsic('toc');

%possible improvements:
%37 point rule from DCUTRI (along with null rules?)
%option to divide into 2 subtriangles like in CUBPACK?

%num evals:  7 16 32 65

%
%***KEYWORDS automatic multidimensional integrator,
%            n-dimensional simplex, general purpose, global adaptive
%***PURPOSE  To calculate an approximation to a vector of integrals
%              over a collection of simplices
%
%               I = I (F ,F ,...,F   ) DS
%                    S  1  2      ncomponents
%
%            where S is a collection ndims-dimensional simplices,
%            and  F = F (X ,X ,...,X  ), K = 1 : 2, ..., ncomponents,
%                  K   K  1  2      ndims
%            and try to satisfy for each component I(K) of I
%              abs( I(K) - INTEGRALS(K) ) < max( ERROR_ABS, ERROR_REL*abs(I(K)) )
%
%  EXAMPLE: n = 4; v = [ zeros(n,1) eye(n) ];
%     [r e] = adsimp(n,v,2,@(x)[x(1)^3;x(3)^4],1000,1e-5); disp([r'; e'])
%

%***DESCRIPTION Computation of integrals over simplical regions.
%            adsimp is a driver for the integration routine SMPSAD,
%            which repeatedly subdivides the region of integration and
%            estimates the integrals and the errors over the subregions
%            with greatest estimated errors until the error request
%            is met or max_F_evals function evaluations have been used.
%
%   ON ENTRY
%
%     ndims     Integer, number of variables. 1 < ndims
%     ncomponents     Integer, number of components of the integral.
%     max_F_evals   Integer, maximum number of F values.
%            NUM_F_EVALS is number F values for each subregion (see RULE),
%            max_F_evals must be >= NUM_SUBREGIONS*NUM_F_EVALS.
%     F      Function  for computing components of the integrand.
%            Input parameter: X, an array that defines the evaluation point.
%            If ncomponents > 1 component, F must output an NFx1 vector.
%     ERROR_ABS Real requested absolute accuracy.
%     ERROR_REL Real requested relative accuracy.
%     RULE    Integer, key to selected local integration rule.
%            RULE = 0 gives the user a (default) degree 7 integration rule.
%            RULE = 1 gives the user a degree 3 integration rule.
%            RULE = 2 gives the user a degree 5 integration rule.
%            RULE = 3 gives the user a degree 7 integration rule.
%            RULE = 4 gives the user a degree 9 integration rule.
%     VRTS Real array of simplex vertices for NUM_SUBREGIONS simplices;
%           the coordinates of vertex J for simplex K must be VRTS(:,J,K).
%
%   ON EXIT
%
%     INTEGRALS  Real array of dimension ncomponents of integral approximations.
%     ABS_ERROR_ESTIMATE  Real array of dimension ncomponents, of absolute accuracy estimates.
%     NUM_F_EVALS_USED Integer, number of FUNSUB calls used by adsimp.
%     FL Integer.
%            FL = 0 for normal exit, when ABS_ERROR_ESTIMATE(K) <=  ERROR_ABS or
%              ABS_ERROR_ESTIMATE(K) <=  abs(INTEGRALS(K))*ERROR_REL with max_F_evals or less
%              function evaluations for all values of K, 1 <= K <= ncomponents.
%            FL = 1 if max_F_evals was too small for adsimp to obtain
%              the required accuracy. In this case adsimp returns
%              values INTEGRALS with estimated absolute accuracies ABS_ERROR_ESTIMATE.
%            FL = 2 if RULE < 0 or RULE > 4,
%            FL = 3 if ndims < 2,
%            FL = 4 if ncomponents < 1,
%            FL = 5 if ERROR_ABS < 0 and ERROR_REL < 0,
%
% if do_time
% adsimp_tic = tic;
% end

[~,~,NUM_SUBREGIONS] = size(VRTS);  %NUM_SUBREGIONS is always 1 since we're only integrating 1 triangle at a time (each triangle has it's own shape function anyway)


switch rule_constants.rule
    case 1
        NUM_F_EVALS = 7;
    case 2
        NUM_F_EVALS = 16;
    case 3
        NUM_F_EVALS = 32;
    case 4
        NUM_F_EVALS = 65;
    otherwise
        NUM_F_EVALS = NaN;
end

[INTEGRALS, ABS_ERROR_ESTIMATE, NUM_F_EVALS_USED, FL] = SMPSAD( ndims, ncomponents,  max_F_evals, ERROR_ABS, ERROR_REL, NUM_F_EVALS, NUM_SUBREGIONS, VRTS, rule_constants, integrand_constants);

%


% if do_time
%     Time = toc(adsimp_tic);
% else
%     Time = [];
% end
