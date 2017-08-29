function printPostSolver(fid, expInfo)
%PRINTPOSTSOLVER   Print commands after solution has been found
%
% Calling sequence:
%   PRINTPOSTSOLVER(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract info from the EXPINFO struct:
allVarNames = expInfo.allVarNames;
allVarString = expInfo.allVarString;
numVars = expInfo.numVars;
indVarNameSpace = expInfo.indVarNameSpace;

% Print commands that will create a plot of the solution obtained:
fprintf(fid, '\n%%%% Plot the solution.\n');
fprintf(fid, 'figure\n');
if ( numVars == 1 ) % Scalar case
    fprintf(fid, 'plot(%s, ''LineWidth'', 2)\n', allVarString);
else
    fprintf(fid, 'plot([%s], ''LineWidth'', 2)\n', allVarString);
end

% Title and xlabel
fprintf(fid, 'title(''IVP solution''), xlabel(''%s'')', indVarNameSpace);

% Deal with ylabel (scalar problem) or legend (systems):
if ( numel(allVarNames) == 1 )
    % Scalar problem:
    fprintf(fid, ', ylabel(''%s'')', allVarNames{:});
else
    % Coupled system. Create a legend.
    leg = '';
    for k = 1:numel(allVarNames)-1
        leg = [leg '''' allVarNames{k} '''' ','];
    end
    leg = [leg '''' allVarNames{k+1} ''''];
    fprintf(fid, ', legend(%s)\n', leg);
end

end
