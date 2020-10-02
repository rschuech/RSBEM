function msg = checkpvpairs(pvpairs,linestyle)
%CHECKPVPAIRS Check length of property value pair inputs
%   MSG = CHECKPVPAIRS(PAIRS) returns an error message if the
%   length of the supplied cell array of property value pairs is
%   incorrect. The message is tailored to accepting LINESPEC as an
%   option convenience input argument.

%   Copyright 1984-2011 The MathWorks, Inc.

msg = '';
npvpairs = length(pvpairs)/2;
if nargin == 1, linestyle = true; end
if (length(pvpairs) == 1) && linestyle
  msg = message('MATLAB:checkpvpairs:UnrecognizedLineStyleParameter');
elseif npvpairs ~= fix(npvpairs)
  msg = message('MATLAB:checkpvpairs:EvenPropValuePairs');
end


