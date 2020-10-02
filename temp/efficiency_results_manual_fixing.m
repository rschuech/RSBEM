
limits = [1 10; 0 1];

for SF1 = 1:0.5:2
    
    % SF1 = 1.5;
    task = 5;
    F_plt = Fs{task};
%     F_plt = F_temp;
%     F_plt = F_fix;
    coord = standardize([SF1 0],limits);  shat = repmat(coord(1),1000,1);  SF2 = linspace(0,1,1000)';
    unstandardized = unstandardize(F_plt.Points,limits);
    y = F_plt([shat SF2]);  figure(829);
    filter = F_plt.Points(:,1) == coord(1);
    SF2_unstd = unstandardize([zeros(size(SF2)) SF2],limits);
    plot(SF2_unstd,y,'-',unstandardized(filter,2),F_plt.Values(filter),'r*'); grid;
    hold on
%     y0 = F_fix([shat SF2]);  filter0 = F_fix.Points(:,1) == coord(1); unstandardized0 = unstandardize(F_fix.Points,limits_uptake);
%      plot(SF2,y0,'-',unstandardized0(filter0,2),F_fix.Values(filter0),'r*'); 
    hold off
    title(['SF1 = ',num2str(SF1)]);
    pause
end