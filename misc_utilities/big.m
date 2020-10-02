function [bigness] = big(var)
%Outputs size in MB of input variable.
bigness = whos('var');
bigness = bigness.bytes / 1E6;