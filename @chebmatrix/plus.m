function C = plus(A, B)
%+         Sum of chebmatrices.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% If an input is not a CHEBMATRIX, need to wrap it in a cell to allow simpler
% code below.
if ( ~isa(A, 'chebmatrix') )
    A = chebmatrix({A});
end
if ( ~isa(B, 'chebmatrix') )
    B = chebmatrix({B});
end

% Set up for addition
[m, n] = size(A);
C = cell(m, n);

% TODO: Surely we want to check that the sizes of A and B are identical.
% Otherwise, the following runs without giving an error:
%   I = operatorBlock.eye([0 1]);
%   A = [ I I ];    B = [I I I];    C = A + B;
% In this case, the C returned as an 1x2 chebmatrix, which can't be what we
% want. Would we still want to allow I + A, where A is e.g. a 1x2 chebmatrix?
% AB, 28/1/14.

% Loop through the blocks of A and B
for i = 1:m
    for j = 1:n
        C{i,j } = A.blocks{i, j} + B.blocks{i, j};
    end
end

% Convert the cell C to a CHEBMATRIX
C = chebmatrix(C);

end
