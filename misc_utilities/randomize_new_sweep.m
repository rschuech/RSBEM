fields = fieldnames(new_sweep);

perm = randperm(length(new_sweep.AR1));

for f = 1:length(fields)
    new_sweep.(fields{f}) = new_sweep.(fields{f})(perm);
end
    