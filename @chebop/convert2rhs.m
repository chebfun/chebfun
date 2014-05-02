function newRHS = convert2rhs(rhs, residual)
%CONVERT2RHS    Convert RHS to a useful format.
%   CHEBOP backslash accepts a variety of syntax for specifying right-hand sides
%   (scalars, CHEBFUNS and CHEBMATRICES). This method takes care of converting
%   the input to a format used internally in chebop.
%
%   This method also takes care of checking that the dimension of the RHS
%   matches the dimensions of the operator. Note that CHEBOP requires the RHS of
%   coupled systems to match the system, even for scalar right-hand sides, e.g.
%       N = chebop(@(x, u, v) [diff(u) + v ; u + diff(v)]);
%       N.bc = @(x, u, v) [u(-1); v(1)];
%       uv = N\0;
%   is not an accepted syntax.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Document inputs.

% Check the size of the residual (the output the dimensions of the CHEBOP).
[numRow, numCol] = size(residual);

% Check whether dimensions match:
if ( ~all(size(rhs) == [numRow, numCol]) )
    if ( all(size(rhs) == [numCol, numRow]) )
        warning('CHEBFUN:CHEBOP:CONVERT2RHS', ...
            'RHS does not match output dimensions of operator.');
        rhs = rhs.';
    else
        error('CHEBFUN:CHEBOP:CONVERT2RHS', ...
            'RHS does not match output dimensions of operator.');
    end
end

% Prepare to convert the RHS to a CHEBMATRIX:
rhsBlocks = cell(numRow, numCol);

% Obtain the blocks of the residual CHEBMATRIX, in order to be able to ensure
% all components of the RHS will be of a type that matches that of the residual.
resBlocks = residual.blocks;

% Need the domain of the residual in order to create the RHS CHEBMATRIX.
dom = residual.domain;

% Convert numerical values in RHS vector into chebmatrix:
for rhsCounter = 1:numRow
    if isa(resBlocks{rhsCounter}, 'chebfun')
        % If corresponding block in the residual is a chebfun, the rhs must also
        % be made to be a chebfun
        rhsBlocks{rhsCounter} = chebfun(rhs(rhsCounter), dom);
    else
        % Otherwise, the entry in the chebmatrix will be a scalar
        rhsBlocks{rhsCounter} = rhs(rhsCounter);
    end
end

% Convert the cell-array to a chebmatrix
newRHS = chebmatrix(rhsBlocks);

end