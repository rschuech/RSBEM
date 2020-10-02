bodies = {...
    [1 0; 1.5 0; 1.5 0.05; 1.5 0.1; 1.5 0.15; 1.5 0.2; 1.5 0.25;...
    1.125 0; 1.125 0.05; 1.125 0.1; 1.125 0.15; 1.125 0.2; 1.125 0.25; 1.125 0.3; 1.125 0.35;]...
    };
    
    [SF1,SF2] = ndgrid(2:0.5:10,[0 0.05]);
    
    bodies{1} = [bodies{1}; [SF1(:) SF2(:)]];
    
    body_tails{1} = [1 0];
    
    bodies{2} = [5 0.2; 5 0.15; 4.5 0.15; 5.5 0.15; 4.5 0.25; 5 0.25; 5.5 0.25;  4.5 0.2; ...
        5.5 0.2; 5.5 0.15; 5 0.15; 6 0.15;  5.5 0.25; 5 0.25; 6 0.25; 6 0.2; 5 0.2;];
    body_tails{2} = [5.5 0.3];
    
    bodies{3} = [4 0.25; 4 0.3;];  body_tails{3} = [4 0.35];
    
    bodies{4} = [3.5 0.3; 3.5 0.35;];  body_tails{4} = [3.5 0.4];
    
    bodies{5} = [3 0.35; 3 0.4;  3 0.45; 3 0.5;];  body_tails{5} = [3 0.55];
    
    bodies{6} = [2.5 0.3; 2.5 0.35;];  body_tails{6} = [2.5 0.25];
    
    bodies{7} = [3 0.75; 3 0.8; 3 0.85; 3 0.9; ...
        2.5 0.75; 2.5 0.775; 2.75 0.775; 2.75 0.8; 2.75 0.825; 2.75 0.85;...
        3.5 0.8; 3.5 0.85; 3.5 0.9; 3.25 0.925; 3.5 0.925; 4 0.95; 4 0.9;];
    body_tails{7} = [3.5 0.7];
    
    bodies{8} = [4 0.4; 4 0.45; 4.5 0.4;]; body_tails{8} = [4.5 0.45];
    
     [SF1,SF2] = ndgrid(7.5:0.5:10, 0.8:0.05:0.95);
     
    bodies{9} = [SF1(:) SF2(:)]; body_tails{9} = [7 0.85];
    
    [SF1,SF2] = ndgrid([6.5:0.5:10],[0.75:-0.05:0.35]);
    
    bodies{10} =  [SF1(:) SF2(:)]; body_tails{10} = [5 0.55];
    
    [SF1,SF2] = ndgrid(7.5:0.5:10,0.25:0.05:0.35);
     bodies{11} =  [SF1(:) SF2(:)]; body_tails{11} = [5 0.55];
    
    bodies{12} = [6 0.4; 6 0.45; 6 0.5; 6 0.55; 6 0.6; 6 0.65; 6 0.7;];
    body_tails{12} = [5 0.55];
    
    bodies{13} = [5.5 0.4; 5.5 0.45; 5.5 0.5; 5.5 0.55; 5.5 0.6; 5.5 0.65;];
    body_tails{13} = [5 0.55];
    
    F_amp = F.temporal.F.amp;  F_amp.ExtrapolationMethod = 'linear';
     F_lambda = F.temporal.F.lambda;  F_lambda.ExtrapolationMethod = 'linear';
      F_nlambda = F.temporal.F.nlambda;  F_nlambda.ExtrapolationMethod = 'linear';
    
    for b = 1:length(bodies)
        tails{b} = [ F_amp(standardize(body_tails{b},limits))  F_lambda(standardize(body_tails{b},limits))  F_nlambda(standardize(body_tails{b},limits))  ];
    end
    
    Bodies = [];  Tails = [];
    for b = [10 11:13 8:-1:1 9]
        
        Bodies = [Bodies; bodies{b}];
        Tails = [Tails; repmat(tails{b},size(bodies{b},1),1)];
    end
    
    [Bodies,ia] = unique(Bodies,'rows','stable');
    Tails = Tails(ia,:);
    
    
    %%
    hold on
    pot = plot(Bodies(:,1),Bodies(:,2),'ko','markersize',6);
    
    
    %%
  
racks = [1  2  3  4];
reps = [1 1 1];
% reps = [1 2 1 2 1 2 1 2];
allbods = [];
for n = 1:length(racks)
    clear sweep
    % n = 4;
    bodies = Bodies(n:length(racks):end,:);
    guesses = Tails(n:length(racks):end,:);
%     inds = randperm(size(bodies,1));
%     bodies = bodies(inds,:);
%     guesses = guesses(inds,:);
sweep.bodies = bodies;  
sweep.amp = guesses(:,1); sweep.lambda = guesses(:,2);  sweep.nlambda = guesses(:,3);

    save(['C:\Hull\SNR_sweep_CFD0',num2str(racks(n)),'.mat'],'sweep');

    
end