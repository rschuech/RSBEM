function [] = whobig()

%displays size of all variables in workspace in MB

%this is a function instead of a script simply to avoid creating useless
%variables (all below variables cleared upon function exit)

temp = evalin('base', 'whos;');

names = {temp.name};
[bytes,inds] = sort([temp.bytes]);

sortnames = names(inds);
sortbigs = bytes / 1E6;  %MB


for bc = 1:length(sortbigs)
    disp([sortnames{bc},'        ',num2str(sortbigs(bc))]);
end

disp(['All Variables:','        ',num2str(sum(sortbigs))]);

