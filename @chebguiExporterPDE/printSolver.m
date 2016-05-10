function printSolver(fid, expInfo)
%PRINTSOLVER   Print commands for solving problems.
%
% Calling sequence:
%   PRINTPOSTSOLVER(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract information from the EXPINFO struct:
allVarString = expInfo.allVarString;
pdeSolver = expInfo.pdeSolver;
indVarName = expInfo.indVarName;
sol0 = expInfo.sol0;

% Add extra whitespace around commas in allVarString
allVarString = strrep(allVarString,',',', ');

% Print commands for solving the problem:
fprintf(fid, '\n%%%% Call %s to solve the problem.\n', pdeSolver);
fprintf(fid, '[%s, %s] = %s(pdefun, %s, %s, bc, opts);\n', indVarName{2}, ...
    allVarString, pdeSolver, indVarName{2}, sol0);

end
