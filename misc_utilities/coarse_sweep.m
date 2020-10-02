
% AR1 = [ repmat( 1.1:0.1:2 , length([0.05:0.05:1]), 1) ];  %refinement in low AR1 region
AR1_0 = 2:1:10;  AR2_0 = 0:0.1:0.9;

AR1 = [ repmat( AR1_0, length(AR2_0), 1) ];  
AR2 = repmat(AR2_0', 1,length(AR1_0));
bods = [ [AR1(:) AR2(:)] ;   ];
[~, obj_signed] = curved_rod_fat_limit(bods(:,1), bods(:,2), 1);
bods = bods(obj_signed > 0, :);
bods = [1 0; bods]; % hard code sphere since fat limit code can't filter out other AR1 = 1, AR2 > 0 impossible shapes

bods = unique( [5 0.5; 1 0; 10 0; 10 0.9; 2 0.6; 3 0.9; 3 0.3; 7 0.3; 4 0.7;  8 0.6; bods(randperm(size(bods,1)),:)], 'stable' , 'rows' );
