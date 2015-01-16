function printSolver(fid, expInfo)
%PRINTSOLVER     Print commands for solving problems
%
% Calling sequence:
%   PRINTPOSTSOLVER(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Print commands for solving the problem:
fprintf(fid,'\n%%%% Solve!\n');
fprintf(fid, ['%% Call solvebvp to solve the problem.\n' ...
    '%% (With the default options, this is equivalent to u = N\\rhs.)\n']);
fprintf(fid, 'u = solvebvp(N, rhs, options);\n');

end