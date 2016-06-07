function h = assignColumns(f, colIdx, g)
%ASSIGNCOLUMNS   Extract columns (or rows) of an array-valued TRIGTECH.
%   G = ASSIGNCOLUMNS(F, COLIDX, G) assigns the columns specified by the row
%   vector COLIDX from the FUN F so that F(:, COLIDX) = G. COLIDX need not be
%   increasing in order or unique, but must contain only integers in the range
%   [1, M] (where F has M columns) and satisfy LENGTH(COLIDX) = SIZE(G, 2) or
%   ISEMPTY(G).
%
% See also EXTRACTCOLUMNS, MAT2CELL.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% g is empty - remove columns:
if ( isempty(g) )
    h = f;
    h.values(:,colIdx) = [];
    h.coeffs(:,colIdx) = [];
    h.isReal(:,colIdx) = [];
    return
end

% Prolong so that f and g have the same length:
if ( length(f) > length(g) )
    g = prolong(g, length(f));
elseif ( length(g) > length(f) )
    f = prolong(f, length(g));
end

% Assign the columns of h.values and h.coeffs:
h = f;
h.values(:, colIdx) = g.values;
h.coeffs(:, colIdx) = g.coeffs;

% Update ishappy:
h.ishappy = f.ishappy && g.ishappy;
h.isReal(colIdx) = g.isReal;

end
