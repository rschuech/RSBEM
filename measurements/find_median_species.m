
temp = vertcat(data_filtered.good.median_unweighted);   SF = vertcat(temp.SF);  med = median(SF)
temp = vertcat(data_filtered.good.mean_unweighted);   SF = vertcat(temp.SF);  avg = median(SF)


[~,inds] = sortrows(SF,[1 2]);


for i = inds'
    disp([data_filtered.good(i).genus,'    ', data_filtered.good(i).species,'    ',num2str(SF(i,:))]);
end