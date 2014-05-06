function out = double2chebmatrix(dVec, residual)
%DOUBLE2CHEBMATRIX  Convert RHS and initial guesses to a useful format.
%   OUT = DOUBLE2CHEBMATRIX(DVEC, RESIDUAL) returns a CHEBMATRIX NEWRHS,
%   constructed from DVEC. Here, DVEC can be a double or a vector of doubles.
%   The input RESIDUAL is a CHEBMATRIX. The dimensions of DVEC must match the
%   dimension of RESIDUAL, that is, if DVEC is a vector, it must have as many
%   components as RESIDUAL. Otherwise, an error is thrown. Note that this method
%   can both be used to convert numerical right-hand sides and numerical initial
%   guesses to the correct format.
%   
%   The kth block of OUT will be an object of the same class as the
%   kth entry of the CHEBMATRIX RESIDUAL.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Developer note:
%   CHEBOP backslash accepts a variety of syntax for specifying right-hand sides
%   (scalars, CHEBFUNS and CHEBMATRICES). This method takes care of converting
%   the input to a format used internally in CHEBOP, that is, a CHEBMATRIX.

% Check the size of the residual (the output the dimensions of the CHEBOP).
[numRow, numCol] = size(residual);


% Prepare to convert the RHS to a CHEBMATRIX:
outBlocks = cell(numRow, numCol);

% Obtain the blocks of the residual CHEBMATRIX, in order to be able to ensure
% all components of the output will be of a type that matches that of the
% residual.
resBlocks = residual.blocks;

% Need the domain of the residual in order to create the RHS CHEBMATRIX.
dom = residual.domain;

% Convert numerical values in RHS vector into a CHEBMATRIX:
for outCounter = 1:numRow
    if isa(resBlocks{outCounter}, 'chebfun')
        % If corresponding block in the residual is a CHEBFUN, the output must
        % also be made to be a CHEBFUN
        outBlocks{outCounter} = chebfun(dVec(outCounter), dom);
    else
        % Otherwise, the entry in the CHEBMATRIX will be a scalar
        outBlocks{outCounter} = dVec(outCounter);
    end
end

% Convert the cell-array to a CHEBMATRIX
out = chebmatrix(outBlocks);

end
