function eps = krondelta(i,j,k)

vec = [i j k];

if isequal(vec,[1 2 3]) ||  isequal(vec,[3 1 2]) ||  isequal(vec,[2 3 1])
    eps = 1;
elseif isequal(vec,[1 3 2]) ||  isequal(vec,[2 1 3]) ||  isequal(vec,[3 2 1])
    eps = -1;
else
    eps = 0;
end