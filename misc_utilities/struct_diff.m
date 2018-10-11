function [diff_fields] = struct_diff(s1,s2)

%outputs the fields that differ between two input structs
%note, probably only works right for comparing different indices of Inputs
%structure

% http://uk.mathworks.com/matlabcentral/fileexchange/22752-compare-structures

[cv_st_1_msgs] = comp_struct(s1,s2,0,1,2*eps,'S1','S2',0);

%[cv_st_1_msgs] = comp_struct(Inputs(1),Inputs(2))

% Declare regular expression to capture field name from comp_struct()
% message output:
regex = '\(.*\)\.(.*)\sand';

for i = 1:numel(cv_st_1_msgs)
    if ~isstr(cv_st_1_msgs{i})  %some numerical values are different, we're screwed
        diff_fields = 'NaN';
        return
    end
end

% Apply (regex) to all comp_struct output strings & return the capture groups as an (N x 1) cell vector of doubly-nested cell vectors containing the token match.
cv_cv_cv_fieldnames = regexp(cv_st_1_msgs,regex,'tokens');


N_fields = numel(cv_cv_cv_fieldnames);
%N_fields = 1;
% Preallocate output:
%cv_fieldnames = cell(N_fields,1);
clear diff_fields
% Iterate through the regexp() output:
ii = 0;
for i_field = 1 : N_fields
if isempty(cv_cv_cv_fieldnames{i_field});
    continue
end
% Assign (i)'th output element:
ii = ii + 1;
diff_fields{ii} = cv_cv_cv_fieldnames{i_field}{1}{1};
end

