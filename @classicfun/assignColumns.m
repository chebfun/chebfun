function f = assignColumns(f, colIdx, g)
%ASSIGNCOLUMNS   Assign columns (or rows) of an array-valued CLASSICFUN.
%   G = ASSIGNCOLUMNS(F, COLIDX) extracts the columns specified by the row
%   vector COLIDX from the CLASSICFUN F so that F(:, COLIDX) = G. COLIDX need
%   not be increasing in order or unique, but must contain only integers in the
%   range [1, M] (where F has M columns) and satisfy LENGTH(COLIDX) = SIZE(G, 2)
%   or ISEMPTY(G).
%
% See also EXTRACTCOLUMNS, MAT2CELL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Assign the columns from the ONEFUN:
if ( isa(g, 'classicfun') )
    f.onefun = assignColumns(f.onefun, colIdx, g.onefun);
else
    f.onefun = assignColumns(f.onefun, colIdx, g);
end

end
