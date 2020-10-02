function [pp, ppder ] = splining(contours)

for i = 1:length(contours)
    if ~isempty(contours(i).x_c)
        t = linspace(0,1,length(contours(i).x_c));
        
        pp(i) = spline(t, [contours(i).x_c; contours(i).y_c]);
        
%         pp(i) = interp1(t,[contours(i).x_c; contours(i).y_c]','linear','pp');
        
%         pp(i) = tpaps(t', [contours(i).x_c;]');
        
        
        ppder(i) = fnder(pp(i));
    else
        pp(i) = NaN;
        ppder(i) = NaN;
    end
    
end