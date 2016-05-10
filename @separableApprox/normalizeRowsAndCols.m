function F = normalizeRowsAndCols(F, p)
%NORMALIZEROWSANDCOLS   Normalize the rows and columns of a SEPARABLEAPPROX.

% Copyright 2016 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

% TODO: Document
% TODO: is this useful?

% Initialise:
numCols = numColumns(F.cols);
colNorms = zeros(1,numCols);
rowNorms = zeros(1,numCols);

% Default to 2-norm:
if ( nargin < 2 )
    p = 2;
end

% Compute norms of each column and row:
for k = 1:numCols
    colNorms(k) = norm(F.cols(:,k), p);
    rowNorms(k) = norm(F.rows(:,k), p);
end

% Scale rows and cols:
F.cols = F.cols./colNorms;
F.rows = F.rows./rowNorms;

% Update pivotValues:
F.pivotValues = F.pivotValues(:).'./colNorms./rowNorms;

end
