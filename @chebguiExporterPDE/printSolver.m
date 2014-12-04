function printSolver(fid, expInfo)
%PRINTSOLVER   Print commands for solving problems.
%
% Calling sequence:
%   PRINTPOSTSOLVER(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract information from the EXPINFO struct:
indVarName = expInfo.indVarName;
sol = expInfo.sol;
sol0 = expInfo.sol0;
deInput = expInfo.deInput;
s = expInfo.s;

% Print commands for solving the problem:
fprintf(fid, '\n%%%% Call pde15s to solve the problem.\n');
fprintf(fid, '[%s, %s] = pde15s(pdefun, %s, %s, bc, opts);\n', indVarName{2}, ...
    sol, indVarName{2}, sol0);

% Conver sol to variable names
if ( numel(deInput) > 1 )
    fprintf(fid, '\n%% Recover variable names.\n');
    for k = 1:numel(s)
        fprintf(fid, '%s = %s(%d,:);\n', s{k}, sol, k);
    end
end

end