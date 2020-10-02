function [INTEGRALS, ABS_ERROR_ESTIMATE, NUM_F_EVALS_USED, FL] = SMPSAD( ndims, ncomponents,  max_F_evals, ERROR_ABS, ERROR_REL, NUM_F_EVALS, NUM_SUBREGIONS, VRTS ,rule_constants,integrand_constants)

%#codegen

%
%***BEGIN PROLOGUE SMPSAD
%***KEYWORDS automatic multidimensional integrator,
%            n-dimensional simplex,
%            general purpose, global adaptive
%
%***LAST MODIFICATION 2012-06
%***PURPOSE  The routine calculates an approximation to a given
%            vector of definite integrals, I, over a hyper-rectangular
%            region hopefully satisfying for each component of I the
%            following claim for accuracy:
%            abs( I(K) - INTEGRALS(K) ) <= max( ERROR_ABS, ERROR_REL*abs(I(K) ) )
%***DESCRIPTION Computation of integrals over hyper-rectangular regions.
%            SMPSAD repeatedly subdivides the regions of integration
%            and estimates the integrals and the errors over the
%            subregions with  greatest estimated errors until the error
%            request is met or maxSUB subregions are stored. The regions
%            are divided into three or four equally sized parts along
%            the direction(s) with greatest absolute fourth difference.
%
%   ON ENTRY
%
%     ndims     Integer, number of variables, ndims > 1.
%     ncomponents     Integer, number of components of the integral.
%     F      Function  for computing components of the integrand.
%            Input parameter: X, an array that defines the evaluation point.
%            If ncomponents > 1 component, F must output an NFx1 vector.
%     max_F_evals Integer.
%            The computations proceed until further subdivision would
%            require more than max_F_evals F values.
%     ERROR_ABS Real, requested absolute accuracy.
%     ERROR_REL Real, requested relative accuracy.
%     RULE    Integer, key to selected local integration rule.
%            RULE = 1 gives the user a degree 3 integration rule.
%            RULE = 2 gives the user a degree 5 integration rule.
%            RULE = 3 gives the user a degree 7 integration rule.
%            RULE = 4 gives the user a degree 9 integration rule.
%     NUM_F_EVALS   Integer, number of FUNSUB calls needed for each subregion.
%     NUM_SUBREGIONS    Integer, number of subregions.
%     VRTS Real array of size (ndims,ndims+1,NUM_SUBREGIONS).
%            Simplex vertices for each subregion; for subregion K vertex
%            J must have components VERTEX(I,J,K), I = 1 : 2, ..., ndims.
%   ON RETURN
%
%     INTEGRALS  Real array of dimension ncomponents.
%            Approximations to all components of the integral.
%     ABS_ERROR_ESTIMATE  Real array of dimension ncomponents.
%            Estimates of absolute accuracies.
%     NUM_F_EVALS_USED Integer, number of new FUNSUB calls used by SMPSAD.
%     FL Integer.
%            FL = 0 for normal exit, when ABS_ERROR_ESTIMATE(K) <=  ERROR_ABS or
%              ABS_ERROR_ESTIMATE(K) <=  abs(INTEGRALS(K))*ERROR_REL, 1 <= K <= ncomponents.
%            FL = 1 if max_F_evals was too small for SMPSAD to obtain the
%            required accuracy. In this case SMPSAD returns values
%            of INTEGRALS with estimated absolute accuracies ABS_ERROR_ESTIMATE.
%
%***end PROLOGUE SMPSAD
%
%   Initialize for rule parameters.
%
NUM_F_EVALS_USED = 0;
DFCOST = 1 + 2*ndims*( ndims + 1 );
%
%     Initialize NUM_F_EVALS_USED, and INTEGRALS and ABS_ERROR_ESTIMATE arrays.
%
INTEGRALS = zeros(ncomponents,1);
ABS_ERROR_ESTIMATE = INTEGRALS;
%
%        Compute weights, generators, PTS.
%

%[G, Weights, PTS] =  SMPRMS( ndims, RULE );  this is always the same so
%run the function outside one time and input G, Weights, PTS into this

FCT = factorial(ndims);
VOL = NaN(1,NUM_SUBREGIONS);
VLS = NaN(ncomponents,NUM_SUBREGIONS);
AES = NaN(ncomponents,NUM_SUBREGIONS);

for K = 1 : NUM_SUBREGIONS  %should always start at 1 for 1 initial triangle but will then increase during adaptive integration
    VOL(K) = abs( det( VRTS(:,1:ndims,K) - VRTS(:,ndims+1,K)*ones(1,ndims) ) ) / FCT;
    
    %
    %     Apply basic rule over each simplex.
    %
    [VLS(:,K), AES(:,K)] = SMPRUL( ndims, VRTS(:,:,K), VOL(K), ncomponents,  rule_constants, integrand_constants);
    INTEGRALS = INTEGRALS + VLS(:,K);
    ABS_ERROR_ESTIMATE = ABS_ERROR_ESTIMATE + AES(:,K);
    NUM_F_EVALS_USED = NUM_F_EVALS_USED + NUM_F_EVALS;
end


% FL = max( ABS_ERROR_ESTIMATE > max(repmat(ERROR_ABS,length(INTEGRALS) / length(ERROR_ABS),1),repmat(ERROR_REL,length(INTEGRALS) / length(ERROR_REL),1).*abs(INTEGRALS)) );
% revised to assume that error tols are same length as integrals
FL = max( ABS_ERROR_ESTIMATE > max(repmat(ERROR_ABS,length(INTEGRALS) / length(ERROR_ABS),1),repmat(ERROR_REL,length(INTEGRALS) / length(ERROR_REL),1).*abs(INTEGRALS)) );

%
%     End initialisation.
%
while ( FL > 0 && NUM_F_EVALS_USED + DFCOST + 4*NUM_F_EVALS <= max_F_evals )
    %
    %     Begin loop while error is too large, and NUM_F_EVALS_USED is not too large.
    %
    %     Adjust INTEGRALS and ABS_ERROR_ESTIMATE.
    %
    [GM, ID] = max(AES);
    if ncomponents > 1
        [GM, ID] = max(GM);
    end
    INTEGRALS = INTEGRALS - VLS(:,ID);
    ABS_ERROR_ESTIMATE = ABS_ERROR_ESTIMATE - AES(:,ID);
    
    %
    %     Determine NEWSBS new subregions.
    %
    
    
    [ VRTS, NEW ] =  SMPDFS( ndims, ID, NUM_SUBREGIONS, VRTS ,integrand_constants);
    
    VI = VOL(ID)/NEW;
    %
    %     Apply basic rule, and add new contributions to INTEGRALS and ABS_ERROR_ESTIMATE.
    %
    K = ID(1);
    VOL(K) = VI;
    [VLS(:,K), AES(:,K)] = SMPRUL( ndims, VRTS(:,:,K), VI, ncomponents,  rule_constants, integrand_constants );
    INTEGRALS = INTEGRALS + VLS(:,K);
    ABS_ERROR_ESTIMATE = ABS_ERROR_ESTIMATE + AES(:,K);
    NUM_F_EVALS_USED = NUM_F_EVALS_USED + NUM_F_EVALS;
    
    VLStemp = NaN(ncomponents,NEW - 1);
    AEStemp = NaN(ncomponents,NEW - 1);
    k = 0;
    for K =   NUM_SUBREGIONS+1:NUM_SUBREGIONS+NEW-1
        k = k + 1;
        [VLStemp(:,k), AEStemp(:,k)] = SMPRUL( ndims, VRTS(:,:,K), VI, ncomponents,  rule_constants, integrand_constants);
        
        INTEGRALS = INTEGRALS + VLStemp(:,k);
        ABS_ERROR_ESTIMATE = ABS_ERROR_ESTIMATE + AEStemp(:,k);
        NUM_F_EVALS_USED = NUM_F_EVALS_USED + NUM_F_EVALS;
    end
    
%     if length(VOL) ~= NUM_SUBREGIONS
%         VOL = NaN;  %screwup check
%     end
    
%     if size(VLS,2) ~= NUM_SUBREGIONS
%         VLS = NaN(size(VLS,1),NUM_SUBREGIONS);  %screwup check
%     end
    
%     if size(AES,2) ~= NUM_SUBREGIONS
%         AES = NaN(size(AES,1),NUM_SUBREGIONS);  %screwup check
%     end
    
    VOL = [VOL repmat(VI,1,NEW - 1)];
    VLS = [VLS VLStemp];
    AES = [AES AEStemp];
    
    
    NUM_F_EVALS_USED = NUM_F_EVALS_USED + DFCOST;
    NUM_SUBREGIONS = NUM_SUBREGIONS + NEW - 1;
    %
    %     Check for error termination.
    %
%     FL = any( ABS_ERROR_ESTIMATE > max(ERROR_ABS,ERROR_REL.*abs(INTEGRALS)) ); % revised to assume that error tols are same length as integrals
   
    FL = max( ABS_ERROR_ESTIMATE > max(repmat(ERROR_ABS,length(INTEGRALS) / length(ERROR_ABS),1),repmat(ERROR_REL,length(INTEGRALS) / length(ERROR_REL),1).*abs(INTEGRALS)) );
    % FL = max( ABS_ERROR_ESTIMATE > max(ERROR_ABS,ERROR_REL*abs(INTEGRALS)) );
    % Coder complains if I don't explicitly expand everything to match sizes
    
    % ERROR_ABS is auto expanded to size of ERROR_REL*abs(INTEGRALS), which
    % is an estimate of requested abs error based on requested rel error.
    % So then the easiest of abs_error or rel_error is taken for each
    % component.  Then these are compared to the actual abs error
    % estimates.  If actual error is smaller, get a 0, else get a 1 for
    % each component.  Finally, taking the max of this vector gives 0 if
    % all errors are OK, otherwise gives a 1 if at least one component has
    % too big of an error
end
%
%     Compute more accurate values of INTEGRALS and ABS_ERROR_ESTIMATE.
%
if NUM_SUBREGIONS > 1
    INTEGRALS = sum(VLS,2);
    ABS_ERROR_ESTIMATE = sum(AES,2);
end
%
%***end SMPSAD