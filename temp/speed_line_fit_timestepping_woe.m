



position = y(:,1:3);

distance = sqrt(   sum(   (position - repmat(position(1,:),size(y,1),1)).^2  ,2)     );


speed = [distance ./ (x - x(1)) ];

inds = 1:1E4:length(distance);

linefits = NaN(length(inds),2);
parfor i = 1:length(inds)
linefits(i,:) = polyfit(x(1:inds(i)),distance(1:inds(i)),1);
end