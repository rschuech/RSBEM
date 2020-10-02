function [ contours] = get_contour(c, x,y,Z)

% x and y must be vectors e.g. unique(X), unique(Y)
% c must be scalar.  need to fix get_contours to handle a vector of c
% values (since each c value can itself have multiple disconnected
% contours)

contours = [];

[temp] = contourcs(x,y,Z,[c c]);

ii = 0;
for i = 1:length(temp)
    if numel(temp(i).X) >= 2  %won't be able to spline with only one data point...
        ii = ii + 1;
        contours(ii).x_c = temp(i).X;
        contours(ii).y_c = temp(i).Y;
    end
    
end
%
% levels = [temp.Level];
%
% for i = 1:length(temp)
%     substruct = temp(levels == c(i));
%     lengths = [substruct.Length];
%     [~,ind] = max(lengths);
%     if ~isempty(ind)
%         contours(i).x_c = substruct(ind).X;
%         contours(i).y_c = substruct(ind).Y;
%     else
%         contours(i).x_c = [];
%         contours(i).y_c = [];
%     end
% end
% 
% 
