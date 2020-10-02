function [ contours] = get_contours(c, x,y,Z)

% x and y must be vectors e.g. unique(X), unique(Y)
if isscalar(c)
    [temp] = contourcs(x,y,Z,[c c]);
else
    [temp] = contourcs(x,y,Z,[c]);
end

levels = [temp.Level];

for i = 1:length(c)
    substruct = temp(levels == c(i));
    lengths = [substruct.Length];
    [~,ind] = max(lengths);
    if ~isempty(ind)
        contours(i).x_c = substruct(ind).X;
        contours(i).y_c = substruct(ind).Y;
    else
        contours(i).x_c = [];
        contours(i).y_c = [];
    end
end


