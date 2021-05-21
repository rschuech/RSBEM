function [] = savedump(filename)


% Get a list of all variables
allvars = evalin('base','whos');
assignin('base','allvars',allvars);
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
assignin('base','tosave',tosave);
% Pass these variable names to save
evalin('base', ['save(', '''', filename ,'''', ', ', 'allvars(tosave).name' ,')' ] );


