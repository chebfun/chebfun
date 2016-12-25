function f = extractColumns(f, colIdx)
%EXTRACTCOLUMNS   Extract columns (or rows) of an array-valued TRIGTECH.
%   G = EXTRACTCOLUMNS(F, COLIDX) extracts the columns specified by the row
%   vector COLIDX from the FUN F so that G = F(:, COLIDX). COLIDX need not be
%   increasing in order or unique, but must contain only integers in the range
%   [1, M], where F has M columns.
%
% See also MAT2CELL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract/re-order the columns from f.values, f.coeffs:
f.values = f.values(:,colIdx);
f.coeffs = f.coeffs(:,colIdx);
f.isReal = f.isReal(colIdx);

end
