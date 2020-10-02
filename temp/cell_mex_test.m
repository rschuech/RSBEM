function [output] = cell_mex_test(input)

celery = cell(input,1);
unit = NaN(10,20);
for i = 1:input
    celery{i} = unit;
end


output = struct('shat',celery);
for i = 1:input
    output(i).shat = rand(10,20);
end