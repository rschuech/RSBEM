
folder = 'C:\Users\rudi\Desktop\RD\body_only_dumps\';

files = dir( [folder , '*AR2_0_forced_dump.mat']);
files = {files.name};

clear AR1 F_coeffs

for f = 1:length(files)
    temp = textscan(   files{f} ,'%s %s %s %f ','Delimiter','_');
    AR1(f) = temp{end};
    load([folder,files{f}],'fcoeffs');
   F_coeffs(f) = fcoeffs;
end

[AR1,inds] = sort(AR1);
F_coeffs = F_coeffs(inds);

save('C:\Users\rudi\Desktop\RD\straight_rod_body_only_friction_coeffs.mat','AR1','F_coeffs');

