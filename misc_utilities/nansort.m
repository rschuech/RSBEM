function [sorted,inds] = nansort(x,mode)

%sorts a vector and outputs sorted and inds that don't include NaNs
if nargin == 2
[sorted, inds] = sort(x,mode);
else
  [sorted, inds] = sort(x);
end

inds = inds(~isnan(sorted));
sorted = sorted(~isnan(sorted));


