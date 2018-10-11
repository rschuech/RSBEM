function [out] = cross2(in1, in2)

% generalized cross prod between in1 and in2 where one is a vector, other
% is a matrix

out = zeros(3,3);

if numel(in1) == 3 %in1 is a vector, cross vector in1 with matrix in2
    for i = 1:3
        for j = 1:3
            
            for k = 1:3
                for l = 1:3
                    out(i,j) = out(i,j) + krondelta(i,k,l) * in1(k) * in2(l,j);
                end
            end
            
        end
    end
    
else % in1 is a matrix, cross matrix in1 with vector in2
    for i = 1:3
        for j = 1:3
            
            for k = 1:3
                for l = 1:3
                    out(i,j) = out(i,j) + in1(i,k) * in2(l) * krondelta(k,l,j);
                end
            end
            
        end
    end
    
end