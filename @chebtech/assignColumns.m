function f = assignColumns(f, colIdx, g)
%ASSIGNCOLUMNS   Extract columns (or rows) of an array-valued CHEBTECH.
%   G = ASSIGNCOLUMNS(F, COLIDX, G) assigns the columns specified by the row
%   vector COLIDX from the FUN F so that F(:, COLIDX) = G. COLIDX need not be
%   increasing in order or unique, but must contain only integers in the range
%   [1, M] (where F has M columns) and satisfy LENGTH(COLIDX) = SIZE(G, 2).
%
% See also EXTRACTCOLUMNS, MAT2CELL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Prolong so that f and g have the same length:
if ( length(f) > length(g) )
    g = prolong(g, length(f));
elseif ( length(g) > length(f) )
    f = prolong(f, length(g));
end

% Assign the columns of f.values and f.coeffs:
f.values(:, colIdx) = g.values;
f.coeffs(:, colIdx) = g.coeffs;

% Update f.vscale:
f.vscale = max(abs(f.values), [], 1);

end
