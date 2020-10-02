function SYMSMS = SMPSMS( ndims, VERTEX, ncomponents,  G, integrand_constants )

%#codegen


%
%***BEGIN PROLOGUE SMPSMS
%***KEYWORDS fully symmetric sum
%***PURPOSE  To compute fully symmetric basic rule sums

%
%***LAST MODIFICATION 2012-06
%***DESCRIPTION SMPSMS computes a fully symmetric sum for a vector
%            of integrand values over a simplex. The sum is taken over
%            all permutations of the generators for the sum.
%
%   ON ENTRY
%
%   ndims       Integer, number of variables.
%   VERTEX  Real array of dimension (1:ndims,1:ndims+1)   2 x 3
%           The vertices of the simplex, one vertex per column.
%   ncomponents      Integer, number of components for the vector integrand.
%   F       Function for computing components of the integrand at X.
%   G       Real Array of dimension (1:ndims+1,1).   3 x 1
%           The generators for the fully symmetric sum.
%
%   ON RETURN
%
%   SYMSMS  Real array of length ncomponents, the values for the fully symmetric
%            sums for each component of the integrand.
%
%***ROUTINES CALLED: Integrand
%
%***end PROLOGUE SMPSMS
%
SYMSMS = zeros(ncomponents,1);
G = -sort(-G);
pr = true;
%
%     Compute integrand value for permutations of G
%

VERTEX = VERTEX(:,:,1);

while pr                             %2 x 3 * 3 x 1
    SYMSMS = SYMSMS + integrand(VERTEX*G, integrand_constants);
    pr = false;
    %
    %     Find next distinct permuation of G and loop back for value.
    %     Permutations are generated in reverse lexicographic order.
    %
    LX = NaN;
    
    for I = 2 : ndims+1
        GI = G(I);
        if G(I-1) > GI
            IX = I - 1;
            for L = 1 : IX/2
                GL = G(L);
                if GL <= GI
                    IX = IX - 1;
                end
                G(L) = G(I-L);
                G(I-L) = GL;
                if G(L) > GI
                    LX = L;
                end
            end
            
            if G(IX) <= GI
                IX = LX;
            end
            G(I) = G(IX);
            G(IX) = GI;
            pr = true;
            break
        end
    end
end
%
%***end SMPSMS
%
%
