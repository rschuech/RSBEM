


V = 1;

%% find max AR2 for given AR1

% AR1 = linspace(1, 7, 2000);  AR1(1) = [];

AR1 = [1.25 1.5 1.75];
AR2 = NaN(size(AR1));


% 1.25  0.35
% 1.75  0.55

for c = 1:length(AR1)
  
    
    objfun = @(AR2) curved_rod_fat_limit(AR1(c),AR2, V);
    
    AR2(c)  = fminbnd(objfun, 0, 1);
    
end

%% find min AR1 for given AR2

AR2 = [0.25 0.5 0.75 0.95];
AR1 = NaN(size(AR2));


for c = 1:length(AR2)
  
    
    objfun = @(AR1) curved_rod_fat_limit(AR1,AR2(c), V);
    
    AR1(c)  = fminbnd(objfun, 1, 10);
    
end