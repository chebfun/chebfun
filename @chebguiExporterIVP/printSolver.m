function printSolver(fid, expInfo)
%PRINTSOLVER     Print commands for solving problems
%
% Calling sequence:
%   PRINTPOSTSOLVER(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Names of variables involved
allVarString = expInfo.allVarString;
numVars = expInfo.numVars;

% Print commands for solving the problem:
fprintf(fid,'\n%%%% Solve!\n');
fprintf(fid, ['%% Call solveivp() to solve the problem.\n' ...
    '%% (With the default options, this is equivalent to u = N\\rhs.)\n']);

if ( numVars == 1) % Scalar case
    fprintf(fid, '%s = solveivp(N, rhs, options);\n', allVarString);
else
    fprintf(fid, '[%s] = solveivp(N, rhs, options);\n', allVarString);
end

end
