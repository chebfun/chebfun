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
sigma = expInfo.sigma;
generalized = expInfo.generalized;

% Print commands for solving the problem:
fprintf(fid, '\n%%%% Solve the eigenvalue problem.\n');

% Need to deal with generalized problems separately.
if ( ~generalized )
    if ( ~isempty(sigma) )
        fprintf(fid, '[V, D] = eigs(N, k, ''%s'', options);\n', sigma);
    else
        fprintf(fid, '[V, D] = eigs(N, k, options);\n');
    end
else
    if ( ~isempty(sigma) )
        fprintf(fid, '[V, D] = eigs(N, B, k, ''%s'', options);\n', sigma);
    else
        fprintf(fid, '[V, D] = eigs(N, B, k, options);\n');
    end
end

end
