function printPostSolver(fid, expInfo)
%PRINTPOSTSOLVER   Print commands after solution has been found.
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
indVarName = expInfo.indVarName;
allVarString = expInfo.allVarString;

% Plotting commands for eigenvalues:
fprintf(fid, '\n%%%% Plot the results.\n');
fprintf(fid, '%% Plot the eigenvalues.\n');
fprintf(fid, 'D = diag(D);\n');
fprintf(fid, 'figure\n');
fprintf(fid, 'plot(real(D), imag(D), ''.'', ''markersize'', 25)\n');
fprintf(fid, 'title(''Eigenvalues''); xlabel(''real''); ylabel(''imag'');\n');

% Plotting commands for eigenmodes:
if ( ischar(allVarNames) || (numel(allVarNames) == 1) )
    fprintf(fid, '\n%% Plot the eigenmodes.\n');
    fprintf(fid, 'figure\n');
    fprintf(fid, 'plot(real(V), ''linewidth'', 2);\n');
    fprintf(fid, 'title(''Eigenmodes''); xlabel(''%s''); ylabel(''%s'');', ...
        indVarName{1}, allVarString);
end

end
