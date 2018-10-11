function [sweep, s] = parameter_sweep(varargin)

%for external function call, just input parameters, which is cell array of cell
%arrays, where first level of organization is parameter type, 2nd level is
%parameter value.  e.g. { {1 2 3}, {4 5 6 7}, {'a' 'b' 'c'} }

% for internal recursive function calls, inputs are 
%(parameters, i, s, sweep)

%output sweep is cell array where row is each parameter
%combo, column is parameter type

if nargin == 1 %initial function call - initialize count variables used during recursion
    
    if isempty(varargin{1})  %actually don't need to sweep anything
        sweep = zeros(1,0);
        s = [];
        return
    end
    
    i = 1;  %counter for parameter type
    s = 1;  %counter for row, i.e. each parameter combo, of output sweep cell array
else
    i = varargin{2};
    s = varargin{3};
    sweep = varargin{4};
end

parameters = varargin{1};

%i is current parameter type at this recursion level

for j = 1:length(parameters{i})  %each value of current parameter type
    
    
    if j > 1 %need to copy previous sweep values for previous parameter types first
        
        for ii = 1:i - 1  %all previous parameter types except the current one
            sweep{s,ii} = sweep{s-1,ii};  %copy last set of values
        end
        
    end
    
    sweep{s,i} = parameters{i}{j};  %insert new value for current parameter type
    
    
    if i < length(parameters)  %if there are still more parameter types, call function again recursively and go to next parameter type
        
        [sweep, s]  = parameter_sweep(parameters, i + 1, s, sweep);  %increment counter for parameter type for next recursion level
        
    else  %if this was the last parameter type, increment row of sweep to move to next parameter combo
        
        s = s + 1; 
        
    end
    
end
