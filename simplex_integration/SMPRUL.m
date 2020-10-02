function [ BASVAL , RGNERR ] = SMPRUL( ndims, VRTS, VOL, ncomponents,  rule_constants, integrand_constants)

%#codegen

VOL = VOL(1);

%
%***BEGIN PROLOGUE SMPRUL
%***KEYWORDS basic numerical integration rule
%***PURPOSE  To compute basic integration rule values.

%***LAST MODIFICATION 2012-06
%***DESCRIPTION SMPRUL computes basic integration rule values for a
%            vector of integrands over a hyper-rectangular region.
%            These are estimates for the integrals. SMPRUL also computes
%            estimates for the errors.
%
%   ON ENTRY
%
%     ndims    Integer, number of variables.
%     VRTS  Real array of size (ndims,ndims+1) of simplex vertices;
%           vertex J must have components VRTS(I,J), I = 1 : 2, ..., ndims.
%     ncomponents    Integer, number of components for the vector integrand.
%     F     Function for computing components of the integrand.
%     VOL   Volume of simplex with VRTS vertices.
%
%   ON EXIT
%
%     BASVAL Real array of length ncomponents, values for the basic rule for
%            each component of the integrand.
%     RGNERR Real array of length ncomponents, error estimates for BASVAL.
%
%***end PROLOGUE SMPRUL
%
RTMN = 1e-1;
SMALL = 1e-12 ;
ERRCOF = 8;
%
%     Compute the rule values.
%
[WTS, RLS] = size(rule_constants.Weights);
RULE = zeros(ncomponents,RLS);

for K = 1 : WTS
    if rule_constants.PTS(K) > 0
               %54x5  1x?      54x1                                                                                                                   %1x5
        RULE = RULE + VOL*SMPSMS( ndims,VRTS,ncomponents, rule_constants.G(:,K), integrand_constants)*rule_constants.Weights(K,:);
    end
end
%
%     Scale integral values and compute the error estimates.
%
BASVAL = zeros(ncomponents,1);
RGNERR = zeros(ncomponents,1);

for I = 1 : ncomponents
    
    BASVAL(I) = RULE(I,1);
    NMBS = abs( BASVAL(I) );
    RT = RTMN;
    NMCP = NaN;
    
    for K = RLS : -2 : 3
        NMRL = max( abs( RULE(I,K) ), abs( RULE(I,K-1) ) );
        
        if ( NMRL > SMALL*NMBS && K < RLS )
            RT = max( NMRL/NMCP, RT );
        end
        
        RGNERR(I) = max( NMRL, RGNERR(I) );
        NMCP = NMRL;
    end
    
    if ( RT < 1 && RLS > 3 )
        RGNERR(I) = RT*NMCP;
    end
    
    RGNERR(I) = max( ERRCOF*RGNERR(I), SMALL*NMBS );
end
%
%***end SMPRUL