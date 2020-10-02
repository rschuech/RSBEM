
U = NaN(size(X));  V = U;  W = U;
parfor i = 1:numel(U)
   
    point = [X(i) Y(i) Z(i)];  %point in full 3D grid
    ind = find(  points(:,1) == point(1) & points(:,2) == point(2) & points(:,3) == point(3) );
    if ~isempty(ind)
        if numel(ind) > 1
            stopafra
        end
%        point
%        ind
    U(i) = field_vel(ind,1);  
    V(i) = field_vel(ind,2);  
    W(i) = field_vel(ind,3);
%    pause
    end
end