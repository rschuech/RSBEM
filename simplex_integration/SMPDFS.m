function  [VRTS, NEW] =  SMPDFS( ndims,  TOP, NUM_SUBREGIONS, VRTS , integrand_constants)

%#codegen


DFNX = NaN;
IE = NaN;
JE = NaN;
LS = 0;
JT = 0;
IT = 0;

TOP = TOP(1);  %coder doesn't know TOP is a scalar

%***BEGIN PROLOGUE SMPDFS
%***PURPOSE  To compute new subregions

%***LAST MODIFICATION 2012-06
%***DESCRIPTION SMPDFS computes fourth differences along each edge
%            direction. It uses these differences to determine a
%            subdivision of the orginal subregion into either three or
%            four new subregions.
%
%   ON ENTRY
%
%   ndims   Integer, number of variables.
%   F    Function for computing the components of the integrand at a point X.
%   TOP    Integer, location in VRTS array for original subregion.
%   NUM_SUBREGIONS Integer, number of subregions in VRTS BEFORE subdivision.
%   VRTS Real array of size (ndims,ndims+1,*), vertices of orginal
%          subregion must be in VRTS(1:ndims,1:ndims+1,TOP).
%
%   ON EXIT
%
%   NEW Integer, number of new subregions (3 or 4).
%   VRTS Real array of size (ndims,1+ndims,*).
%          The vertices of the of new subegions will be at locations
%          TOP, NUM_SUBREGIONS+1 : ..., NUM_SUBREGIONS+NEW-1.
%
CUTTF = 2;
CUTTB = 8;
%
%       Compute the differences.
%
IS = 1;
JS = 2;
DFMX = 0;
EMX = 0;
V = VRTS(:,:,TOP);
CN = mean(V,2);
FC = integrand(CN, integrand_constants);
DFMD = sum( abs(FC) );
FRTHDF = zeros(ndims,ndims+1);

for I = 1 : ndims
    for J = I+1 : ndims+1
        H = 2*( V(:,I)-V(:,J) ) / ( 5*(ndims+1) );
        
        EWD = sum(abs(H));
        if EWD >= EMX
            IE = I;
            JE = J; 
            EMX = EWD;
        end
        
        DFR = sum( abs( integrand(CN-2*H, integrand_constants)+integrand(CN+2*H, integrand_constants) + 6*FC - 4*(integrand(CN-H, integrand_constants)+integrand(CN+H, integrand_constants)) ) );
        DFR = DFR(1); % appease Coder (pretty sure should always be scalar?)
        if DFMD + DFR/8 == DFMD
            DFR = 0;
        end
        DFR = DFR*EWD;
        if DFR >= DFMX
            IT = IS;
            JT = JS;
            DFNX = DFMX;
            IS = I;
            JS = J;
            DFMX = DFR; 
        elseif DFR >= DFNX
            IT = I;
            JT = J;
            DFNX = DFR;
        end
        FRTHDF(I,J) = DFR;
    end
end
%
%     Determine whether to compute three or four new subregions.
%
if DFNX > DFMX/CUTTF
    NEW = 4;
else
    NEW = 3;
    if DFMX == 0
        IS = IE;
        JS = JE;
    else
        DFSMX = 0;
        for L = 1 : ndims+1
            if L ~= IS && L ~= JS
                IT = min( [L IS JS] );
                JT = max( [L IS JS] );
                LT = IS + JS + L - IT - JT;
                DFR =  FRTHDF(IT,LT) + FRTHDF(LT,JT);
                if DFR >= DFSMX
                    DFSMX = DFR;
                    LS = L;
                end
            end
        end
        DIFIL = FRTHDF( min(IS,LS), max(IS,LS) );
        DIFLJ = FRTHDF( min(JS,LS), max(JS,LS) );
        DFNX = DIFIL + DIFLJ - min( DIFIL,DIFLJ );
        if DFMX/CUTTB < DFNX && DIFIL > DIFLJ
            IT = IS;
            IS = JS;
            JS = IT;
        end
    end
end
%
%     Copy vertices and volume for TOP to new subregions

% if size(VRTS,3) ~= NUM_SUBREGIONS
%     VRTS = NaN(2,3);   %hopefully make code error / fail if this ever happens (not sure it ever would)
% end
VRTS = cat(3,VRTS,repmat(V,1,1,NEW-1));

VTI = V(:,IS);
VTJ = V(:,JS);
if ( NEW == 4 ) %     Compute four new subregions.
    VRTS(:,JS,TOP)   = ( VTI + VTJ )/2;
    VRTS(:,IS,NUM_SUBREGIONS+1) = VTI;
    VRTS(:,JS,NUM_SUBREGIONS+1) = ( VTI + VTJ )/2;
    VRTS(:,IS,NUM_SUBREGIONS+2) = ( VTI + VTJ )/2;
    VRTS(:,JS,NUM_SUBREGIONS+2) = VTJ;
    VRTS(:,IS,NUM_SUBREGIONS+3) = ( VTI + VTJ )/2;
    VRTS(:,JS,NUM_SUBREGIONS+3) = VTJ;
    VTI = VRTS(:,IT,TOP);
    VTJ = VRTS(:,JT,TOP);
    VRTS(:,JT,TOP)   = ( VTI + VTJ )/2;
    VRTS(:,IT,NUM_SUBREGIONS+1) = ( VTI + VTJ )/2;
    VRTS(:,JT,NUM_SUBREGIONS+1) = VTJ;
    VTI = VRTS(:,IT,NUM_SUBREGIONS+2);
    VTJ = VRTS(:,JT,NUM_SUBREGIONS+2);
    VRTS(:,JT,NUM_SUBREGIONS+2) = ( VTI + VTJ )/2;
    VRTS(:,IT,NUM_SUBREGIONS+3) = ( VTI + VTJ )/2;
    VRTS(:,JT,NUM_SUBREGIONS+3) = VTJ;
else %     Compute three new subregions.
    VRTS(:,JS,TOP)   = ( 2*VTI + VTJ )/3;
    VRTS(:,IS,NUM_SUBREGIONS+1) = ( 2*VTI + VTJ )/3;
    if ( DFMX/CUTTF < DFNX )
        VRTS(:,JS,NUM_SUBREGIONS+1) = VTJ;
        VRTS(:,IS,NUM_SUBREGIONS+2) = ( 2*VTI + VTJ )/3;
        VRTS(:,JS,NUM_SUBREGIONS+2) = VTJ;
        VTJ = VRTS(:,JS,NUM_SUBREGIONS+1);
        VTL = VRTS(:,LS,NUM_SUBREGIONS+1);
        VRTS(:,LS,NUM_SUBREGIONS+1) = ( VTJ + VTL )/2;
        VRTS(:,JS,NUM_SUBREGIONS+2) = ( VTJ + VTL )/2;
        VRTS(:,LS,NUM_SUBREGIONS+2) = VTL;
    else
        VRTS(:,JS,NUM_SUBREGIONS+1) = ( VTI + 2*VTJ )/3;
        VRTS(:,IS,NUM_SUBREGIONS+2) = ( VTI + 2*VTJ )/3;
        VRTS(:,JS,NUM_SUBREGIONS+2) = VTJ;
    end
end
%
%***end SMPDFS
%
