
SI_file = 'C:\allsync\all papers\curved rod ms\Table S2.xlsx';

SI_data = [];
for d = 1:length(data_filtered.good)
    SI_data{d,1} = data_filtered.good(d).genus;
    SI_data{d,2} = data_filtered.good(d).species;
    SI_data{d,3} = data_filtered.good(d).n_individuals;
    SI_data{d,4} = data_filtered.good(d).median_unweighted.SF(1);  SI_data{d,5} = data_filtered.good(d).median_unweighted.SF(2);
    SI_data{d,6} = data_filtered.good(d).median_unweighted.sph_dia;
end

    
xlswrite(SI_file,SI_data,'data','A3');

