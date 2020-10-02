

AR1 = [1:0.125:2];
AR2 = 0.025:0.025:1;
[AR1, AR2] = ndgrid(AR1,AR2);
AR1AR2 = [AR1(:) AR2(:)];

AR1 = 2:0.25:5;
AR2 = 0.3:0.025:1;
[AR1, AR2] = ndgrid(AR1,AR2);
AR1AR2 = [  AR1AR2; [AR1(:) AR2(:)] ];

AR1AR2 = setdiff(AR1AR2, Best.eff.body, 'rows');

[obj, obj_signed] = curved_rod_fat_limit(AR1AR2(:,1), AR1AR2(:,2), 1);
AR1AR2(obj_signed < 0,:) = [];

AR1AR2(AR1AR2(:,1) == 1 & AR1AR2(:,2) > 0,:) = [];

AR1AR2 = roundn(AR1AR2,-8);