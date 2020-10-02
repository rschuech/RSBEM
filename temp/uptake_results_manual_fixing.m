dump = load('C:\Hull\uptake shat.mat');

% limits = [1 10; 0 0.985];
limits = [1 10; 0 1]; % now it's the same as everything else (was above definition in original stored interpolant)

standardized = standardize(dump.AR1_AR2,limits);

F_fix = scatteredInterpolant(standardized(:,1),standardized(:,2),dump.flux,'linear','none');

%% plot uptake vs SF2

% limits_uptake = [1 10; 0 0.985];
limits_uptake = [1 10; 0 1];

for SF1 = 1:0.1:2
    
    % SF1 = 1.5;
    task = 1;
%     F_plt = Fs{task};
    F_plt = F_temp;
%     F_plt = F_fix;
    coord = standardize([SF1 0],limits_uptake);  shat = repmat(coord(1),3000,1);  SF2 = linspace(0,1,3000)';
    unstandardized = unstandardize(F_plt.Points,limits_uptake);
    y = F_plt([shat SF2]);  figure(829);
    filter = F_plt.Points(:,1) == coord(1);
    plot(SF2,y,'-',unstandardized(filter,2),F_plt.Values(filter),'r*'); grid;
    hold on
    y0 = F_fix([shat SF2]);  filter0 = F_fix.Points(:,1) == coord(1); unstandardized0 = unstandardize(F_fix.Points,limits_uptake);
     plot(SF2,y0,'-',unstandardized0(filter0,2),F_fix.Values(filter0),'r*'); 
    hold off
    title(['SF1 = ',num2str(SF1)]);
    pause
end

%% enforce that uptake decreases with SF2
options = slmset('decreasing','on','knots',-1);
% F_temp = Fs{1};
% F_temp = F_fix;
bads = [1.1 0.1; 1.2 0.3; 1.2 0.35; 1.2 0.375; 1.3 0.375; ];
bads_std = standardize(bads,limits_uptake);

SF1_stds = unique(F_temp.Points(:,1));
c = 0;  new_values = [];
for SF1_std = SF1_stds'
    c = c+1;
    inds = ismember(F_temp.Points(:,1),SF1_std);
    SF2 = F_temp.Points(inds,2);  metric = F_temp.Values(inds); 
    if sum(inds) >= 3
    
    [SF2_temp,inds2] = sort(SF2);  metric = metric(inds2);
    pts_temp = [repmat(SF1_std,length(SF2_temp),1)  SF2_temp];
    [~,ia] = setdiff(pts_temp,bads_std,'rows');
    slm = slmengine(SF2_temp(ia),metric(ia),options);
    
    F_temp.Values( find(inds) ) = slmeval( SF2 , slm ); %apparently can't use logical indexing with a scatteredInterpolant object
%     new_values = [new_values; [SF2  slmeval(SF2,slm)]  ];
    else
%         new_values = [new_values; [SF2 metric] ];
    end
    
end

%% enforce that uptake increases with SF1
options = slmset('increasing','on','knots',-1);
% F_temp = Fs{1};
F_temp = F_fix;
bads = [1.1 0.1; 1.2 0.3; 1.2 0.35; 1.2 0.375; 1.3 0.375; ];
bads_std = standardize(bads,limits_uptake);

SF2_stds = unique(F_temp.Points(:,2));
c = 0;  new_values = [];
for SF2_std = SF2_stds'
    c = c+1;
    inds = ismember(F_temp.Points(:,2),SF2_std);
    SF1 = F_temp.Points(inds,1);  metric = F_temp.Values(inds); 
    if sum(inds) >= 3
    
    [SF1_temp,inds2] = sort(SF1);  metric = metric(inds2);
    pts_temp = [SF1_temp   repmat(SF2_std,length(SF1_temp),1)  ];
    [~,ia] = setdiff(pts_temp,bads_std,'rows');
    slm = slmengine(SF1_temp(ia),metric(ia),options);
    
    F_temp.Values( find(inds) ) = slmeval( SF1 , slm ); %apparently can't use logical indexing with a scatteredInterpolant object
%     new_values = [new_values; [SF2  slmeval(SF2,slm)]  ];
    else
%         new_values = [new_values; [SF2 metric] ];
    end
    
end
    
%% plot uptake vs SF1
for SF2 = 0:0.01:0.5
    
    % SF1 = 1.5;
    task = 1;
%     F_plt = Fs{task};
    F_plt = F_temp;
    coord = standardize([5 SF2],limits_uptake);  shat = repmat(coord(2),1000,1);  SF1_std = linspace(0,0.5,1000)';
    y = F_plt([SF1_std shat]);  figure(832);
    filter = F_plt.Points(:,2) == coord(2);
    plot(SF1_std,y,'-',F_plt.Points(filter,1),F_plt.Values(filter),'r*'); grid;
    xlim([0 0.05])
    pause
end




   