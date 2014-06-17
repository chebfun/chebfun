function F = normalizePivots(F)
%NORMALIZEPIVOTS   Scale rows and cols of a CHEBFUN2 so that all pivots are 1.
%
% Additionally, the norm of the kth row and column will be the same.

% Copyright 2014 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

% TODO: Document
% TODO: is this useful?

F = normalizeRowsAndCols(F);

d = F.pivotValues(:).';
s = sign(d);
sqrtp = sqrt(abs(d));
F.cols = F.cols./(s.*sqrtp);
F.rows = F.rows./sqrtp;
F.pivotValues = ones(1, length(d));

end
