rack = 4;

load('C:\Users\rudi\Desktop\RD\repeats.mat')

for i = 1:length(Repeats)
    if Repeats(i).rack == rack
        movefile(['C:\Users\rudi\Desktop\RD\swept_dumps\',Repeats(i).name,'_*'],'C:\Users\rudi\Desktop\RD\repeats\');
    end
  
end