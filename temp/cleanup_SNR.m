load('C:\Users\rudi\Desktop\RD\Results\Results_Temporal_SN_tryagain_1_synced2.mat')

newbods = unique( [  [Results.AR1]'  [Results.AR2]'  ] , 'rows');

load('C:\Users\rudi\Desktop\RD\Results\Results_Temporal_SN_radial_guesses_1_synced2.mat')

oldbods =   [  [Results.AR1]'  [Results.AR2]'  ] ;

[a,b] = ismember(oldbods,newbods,'rows');

Results(a) = [];
save('C:\Users\rudi\Desktop\RD\Results\Results_Temporal_SN_radial_guesses_1_synced2_cleanedagain.mat','Results')
