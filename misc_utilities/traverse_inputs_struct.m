function [fieldpaths, values] = traverse_inputs_struct(varargin)

%external function call inputs a struct with nested structs

%outputs fieldpaths, a cell array of the "paths" through the nested struct
%to get to the data-containing fields, and values, a cell array of the
%field values

%internal recursive function calls use additional inputs, hence use of
%varargin

if nargin == 1 %first function call, initialize other variables used in recursion
    fieldpath = {};
    fieldpaths = [];
    values = [];
else
    fieldpath = varargin{2};
    fieldpaths = varargin{3};
    values = varargin{4};
end


S = varargin{1};

fields = fieldnames(S);

for f = 1:length(fields)
    currentfield = fields{f};
    
    fieldpath = [fieldpath,currentfield ];
    
    if isstruct(S.(currentfield))  %a nested structure - recurse into it
        
        
        [fieldpaths,values] = traverse_inputs_struct(S.(currentfield),fieldpath,fieldpaths,values);
        
        
    else  %must just be a single value / single array of values to use for all sweep iterations
        
        value = S.(currentfield);
        
        
        fieldpaths{end+1} = fieldpath; %append to cell array of full fieldpaths
        values{end+1} = value;  %append to cell array of parameter values
        
    end
    
    fieldpath(end) = [];  %I don't remember why this is here but apparently it's pretty important
    
end