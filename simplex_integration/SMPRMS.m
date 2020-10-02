function [G, Weights, PTS] =  SMPRMS( ndims, RULE )
%
%***BEGIN PROLOGUE SMPRMS
%***KEYWORDS basic integration rule, degree 2*RULE+1
%***PURPOSE  To initialize a degree 2*RULE+1 basic rule and null rules.

%
%***LAST MODIFICATION 2012-06
%***DESCRIPTION  SMPRMS initializes a degree 2*RULE+1 rule, and
%                and max(2*RULE,2) lower degree null rules.
%
%   ON ENTRY
%
%   ndims    Integer, number of variables.
%   RULE    Integer, <= 4 and >= 0, rule parameter.
%          If RULE > 0 a degree 2*RULE+1 rule is initialized.
%
%   ON Exit
%
%   Weights      Real array of weights for the basic and null rules.
%          Weights(1,1),...,Weights(WTS,1) are weights for the basic rule.
%          Weights(I,1),...,Weights(WTS,I) for I > 1 are null rule weights.
%   G      Real array of fully symmetric sum generators for the rules.
%          G(1,J), ..., G(ndims+1,J) are the generators for the
%          points associated with the Jth weights.
%   PTS    Integer array numbers of integrand
%          values needed for generator J.
%
%***REFERENCES
%
%  Axel Grundmann and H. M. Moller
%  "Invariant Integration Formulas for the n-Simplex by Combinatorial
%    Methods", SIAM J Numer. Anal. 15(1978), 282--290,
% and
%  A. H. Stroud
%  "A Fifth Degree Integration Formula for the n-Simplex
%  SIAM J Numer. Anal. 6(1969), 90--98,
% and
%  I. P. Mysovskikh
%  "On a cubature formula for the simplex"
%  Vopros. Vycisl. i Prikl. Mat., Tashkent 51(1978), 74--90.
%
%
%***ROUTINES CALLED NONE
%***end PROLOGUE SMPRMS
%
%
%     Initialize RLS and GMS.
%
switch RULE
    case 1
        RLS = 3;
        GMS =  2;
        WTS =  3;
    case 2
        RLS = 5;
        GMS =  4;
        WTS =  6;
    case 3
        RLS = 7;
        GMS =  7;
        WTS = 11;
    case 4
        RLS = 7;
        
        if ndims == 2
            GMS = 11;
            WTS = 20;
        else
            GMS = 12;
            WTS = 21;
        end
end

Weights = zeros(WTS,RLS);  %WTS is number weights (so number of points?) and RLS is number of null rules - 1 ?
PTS = zeros(WTS,1);  %each quad point?
G = zeros(ndims+1,WTS);
%
%     Compute generator, PTS and weight values for all rules.
%
NP = ndims + 1;
N2 = NP*( ndims + 2 );
N4 = N2*( ndims + 3 )*( ndims + 4 );
N6 = N4*( ndims + 5 )*( ndims + 6 );
N8 = N6*( ndims + 7 )*( ndims + 8 );

G(:,1) = 1/NP;
PTS(1) = 1;
R1 = ( ndims + 4 - sqrt(15) )/( ndims*ndims + 8*ndims + 1 );
S1 = 1 - ndims*R1;
L1 = S1 - R1;
G(1,GMS+1) = S1;
G(2:NP,GMS+1) = R1;
PTS(GMS+1) = ndims + 1;
IW = RLS;

if ( RULE < 4 )
    %
    %        Compute weights for special degree 1 rule.
    %
    Weights(1,IW) = 1;
    IW = IW - 1;
    Weights(GMS+1,IW) = 1/NP;
    IW = IW - 1;
end
%
%     Compute weights, generators and PTS for degree 3 rule.
%
G(1,2) = 3/( ndims + 3 );
G(2:NP,2) = 1/( ndims + 3 );
PTS(2) = NP;
Weights(2,IW) = ( ndims + 3 )^3/( 4*N2*( ndims + 3 ) );

if ( RULE > 1 )
    IW = IW - 1;
    %
    %        Compute weights, generators and PTS for degree 3 and 5 rules.
    %
    if ( ndims == 2 )
        %
        %           Special degree 3 rule.
        %
        L2 = .62054648267200632589046034361711;
        L1 = -sqrt( 1/2 - L2^2 );
        R1 = ( 1 - L1 )/3;
        S1 = 1 - 2*R1;
        G(1,GMS+1) = S1;
        G(2:NP,GMS+1) = R1;
        PTS(GMS+1) = 3;
        Weights(GMS+1,IW) = 1/6;
        R2 = ( 1 - L2 )/3;
        S2 = 1 - 2*R2;
        G(1,GMS+2) = S2;
        G(2:NP,GMS+2) = R2;
        PTS(GMS+2) = 3;
        Weights(GMS+2,IW) = 1/6;
    else
        %
        %           Degree 3 rule using Stroud points.
        %
        R2 = ( ndims + 4 + sqrt(15) )/( ndims*ndims + 8*ndims + 1 );
        S2 = 1 - ndims*R2;
        L2 = S2 - R2;
        G(1,GMS+2) = S2;
        G(2:NP,GMS+2) = R2;
        PTS(GMS+2) = NP;
        Weights(GMS+2,IW) = ( 2/(ndims+3) - L1 )/( N2*(L2-L1)*L2^2 );
        Weights(GMS+1,IW) = ( 2/(ndims+3) - L2 )/( N2*(L1-L2)*L1^2 );
    end,
    IW = IW - 1;
    %
    %        Grundmann-Moller degree 5 rule.
    %
    G(1,3) = 5/( ndims + 5 );
    G(2:NP,3) = 1/( ndims + 5 );
    PTS(3) = NP;
    G(1,4) = 3/( ndims + 5 );
    G(2,4) = 3/( ndims + 5 );
    G(3:NP,4) = 1/( ndims + 5 );
    PTS(4) = NP*ndims/2;
    Weights(2,IW) = -( ndims + 3 )^5/( 16*N4 );
    Weights(3:4,IW) =  ( ndims + 5 )^5/( 16*N4*( ndims + 5 ) );
end

if ( RULE > 2 )
    IW = IW - 1;
    %
    %        Compute weights, generators and PTS for degree 5 and 7 rules.
    %
    %
    %        Stroud degree 5 rule.
    %
    U1 = ( ndims + 7 + 2*sqrt(15) )/( ndims*ndims + 14*ndims - 11 );
    V1 = ( 1 - ( ndims - 1 )*U1 )/2;
    D1 = V1 - U1;
    G(1,GMS+3) = V1;
    G(2,GMS+3) = V1;
    G(3:NP,GMS+3) = U1;
    PTS(GMS+3) = ( ( ndims + 1 )*ndims )/2;
    U2 = ( ndims + 7 - 2*sqrt(15) )/( ndims*ndims + 14*ndims - 11 );
    V2 = ( 1 - ( ndims - 1 )*U2 )/2;
    D2 = V2 - U2;
    G(1,GMS+4) = V2;
    G(2,GMS+4) = V2;
    G(3:NP,GMS+4) = U2;
    PTS(GMS+4) = ( ( ndims + 1 )*ndims )/2;
    if ( ndims == 2 )
        Weights(GMS+3,IW) = ( 155 - sqrt(15) )/1200;
        Weights(GMS+4,IW) = ( 155 + sqrt(15) )/1200;
        Weights(1,IW) = 1 - 3*( Weights(GMS+3,IW) + Weights(GMS+4,IW) ) ;
    elseif ( ndims == 3 )
        Weights(GMS+1,IW) = ( 2665 + 14*sqrt(15) )/37800;
        Weights(GMS+2,IW) = ( 2665 - 14*sqrt(15) )/37800;
        Weights(GMS+3,IW) = 2*15/567;
        PTS(GMS+4) = 0;
    else
        Weights(GMS+1,IW) = ( 2*(27-ndims)/(ndims+5)-L2*(13-ndims) )/( L1^4*(L1-L2)*N4 );
        Weights(GMS+2,IW) = ( 2*(27-ndims)/(ndims+5)-L1*(13-ndims) )/( L2^4*(L2-L1)*N4 );
        Weights(GMS+3,IW)=( 2/( ndims + 5 ) - D2 )/( N4*( D1 - D2 )*D1^4 );
        Weights(GMS+4,IW)=( 2/( ndims + 5 ) - D1 )/( N4*( D2 - D1 )*D2^4 );
    end,
    IW = IW - 1;
    %
    %        Grundmann-Moller degree 7 rule.
    %
    G(1,5) = 7/( ndims + 7 );
    G(2:NP,5) = 1/( ndims + 7 );
    PTS(5) = NP;
    G(1,6) = 5/( ndims + 7 );
    G(2,6) = 3/( ndims + 7 );
    G(3:NP,6) = 1/( ndims + 7 );
    PTS(6) = NP*ndims;
    G(1:3,7) = 3/( ndims + 7 );
    G(4:NP,7) = 1/( ndims + 7 );
    PTS(7) = NP*ndims*( ndims - 1 )/6;
    Weights(2  ,IW) =  ( ndims + 3 )^7/( 2*64*N4*( ndims + 5 ) );
    Weights(3:4,IW) = -( ndims + 5 )^7/(   64*N6 );
    Weights(5:7,IW) =  ( ndims + 7 )^7/(   64*N6*( ndims + 7 ) );
end
if ( RULE == 4 )
    IW = IW - 1;
    %
    %        Compute weights, generators and PTS for degree 7, 9 rules.
    %
    %        Mysovskikh degree 7 rule.
    %
    SG = 1/( 23328*N6 );
    U5 = -6^3*SG*( 52212 - ndims*( 6353 + ndims*( 1934 - ndims*27 ) ) );
    U6 =  6^4*SG*(  7884 - ndims*( 1541 - ndims*9 ) );
    U7 = -6^5*SG*(  8292 - ndims*( 1139 - ndims*3 ) )/( ndims + 7 );
    P0 = -144*( 142528 + ndims*( 23073 - ndims*115 ) );
    P1 = -12*( 6690556 + ndims*( 2641189 + ndims*( 245378 - ndims*1495 ) ) );
    P2 = -16*(6503401 + ndims*(4020794+ndims*(787281+ndims*(47323-ndims*385))));
    P3 = -( 6386660 + ndims*(4411997+ndims*(951821+ndims*(61659-ndims*665))) )*( ndims + 7 );
    A = P2/( 3*P3 );
    P = A*( P1/P2 - A );
    Q = A*( 2*A*A - P1/P3 ) + P0/P3;
    R = sqrt(-P^3);
    TH = acos(-Q/(2*R))/3;
    R = 2*R^(1/3);
    TP = 2*pi/3;
    A1 = -A + R*cos(TH);
    A2 = -A + R*cos(TH+2*TP);
    A3 = -A + R*cos(TH+TP);
    G(1,GMS+5) = ( 1 - ndims*A1 )/NP;
    G(2:NP,GMS+5) = ( 1 + A1 )/NP;
    PTS(GMS+5) = NP;
    G(1,GMS+6) = ( 1 - ndims*A2 )/NP;
    G(2:NP,GMS+6) = ( 1 + A2 )/NP;
    PTS(GMS+6) = NP;
    G(1,GMS+7) = ( 1 - ndims*A3 )/NP;
    G(2:NP,GMS+7) = ( 1 + A3 )/NP;
    PTS(GMS+7) = NP;
    Weights(GMS+5,IW) = ( U7-(A2+A3)*U6+A2*A3*U5 )/( A1^2-(A2+A3)*A1+A2*A3 )/A1^5;
    Weights(GMS+6,IW) = ( U7-(A1+A3)*U6+A1*A3*U5 )/( A2^2-(A1+A3)*A2+A1*A3 )/A2^5;
    Weights(GMS+7,IW) = ( U7-(A2+A1)*U6+A2*A1*U5 )/( A3^2-(A2+A1)*A3+A2*A1 )/A3^5;
    G(1,GMS+8) = 4/( ndims + 7 );
    G(2,GMS+8) = 4/( ndims + 7 );
    G(3:NP,GMS+8) = 1/( ndims + 7 );
    PTS(GMS+8) = NP*ndims/2;
    Weights(GMS+8,IW) = 10*(ndims+7)^6/( 729*N6 );
    G(1,GMS+9) = 11/( ndims + 7 )/2;
    G(2,GMS+9) =  5/( ndims + 7 )/2;
    G(3:NP,GMS+9) = 1/( ndims + 7 );
    PTS(GMS+9) = NP*ndims;
    Weights(GMS+9,IW) = 64*(ndims+7)^6/( 6561*N6 );
    Weights(4,IW) = Weights(4,IW+1);
    Weights(7,IW) = Weights(7,IW+1);
    IW = IW - 1;
    %
    %        Grundmann-Moller degree 9 rule.
    %
    G(1,8) = 9/( ndims + 9 );
    G(2:NP, 8) = 1/( ndims + 9 );
    PTS(8) = NP;
    G(1,9) = 7/( ndims + 9 );
    G(2,9) = 3/( ndims + 9 );
    G(3:NP, 9) = 1/( ndims + 9 );
    PTS(9) = NP*ndims;
    G(1:2,10) = 5/( ndims + 9 );
    G(3:NP,10) = 1/( ndims + 9 );
    PTS(10) = NP*ndims/2;
    G(1,11) = 5/( ndims + 9 );
    G(2:3,11) = 3/( ndims + 9 );
    G(4:NP,11) = 1/( ndims + 9 );
    PTS(11) = NP*ndims*( ndims - 1 )/2;
    Weights(2   ,IW) = -( ndims + 3 )^9/( 6*256*N6 );
    Weights(3:4 ,IW) =  ( ndims + 5 )^9/( 2*256*N6*( ndims + 7 ) );
    Weights(5:7 ,IW) = -( ndims + 7 )^9/(   256*N8 );
    Weights(8:11,IW) =  ( ndims + 9 )^9/(   256*N8*( ndims + 9 ) );
    
    if ( ndims > 2 )
        G(1:4,12) = 3/( ndims + 9 );
        G(5:NP,12) = 1/( ndims + 9 );
        PTS(12) = NP*ndims*( ndims - 1 )*( ndims - 2 )/24;
        Weights(12,IW) = Weights(8,IW);
    end
    
end
%
%     Compute unnormalized weights.
%
Weights(1,:) = 1 - PTS(2:WTS)'*Weights(2:WTS,:);
NB = PTS'*Weights(:,1).^2;
Weights(:,2:RLS) = Weights(:,2:RLS) - Weights(:,1)*ones(1,RLS-1);
%
%        Orthogonalize and normalize null rules.
%
for K = 2 : RLS
    Weights(:,K) = Weights(:,K) - Weights(:,2:K-1)*Weights(:,2:K-1)'*( PTS.*Weights(:,K) )/NB;
    Weights(:,K) = Weights(:,K)*sqrt( NB/sum(PTS.*Weights(:,K).^2) );
end
%
%***end SMPRMS
%