function f = assignColumns(f, colIdx, g)
%ASSIGNCOLUMNS   Assign columns (or rows) of an array-valued CHEBFUN.
%   F = ASSIGNCOLUMNS(F, COLIDX, G), if F is a column CHEBFUN, assigns the
%   columns of the CHEBFUN G (or rows, if G is a row CHEBFUN) to the columns
%   specified by the vector COLIDX so that F(:,COLIDX) = G. COLIDX need not be
%   increasing in order or unique but must contain only integers in the range
%   [1, M] (where F has M columns) and have a length equal to the number of
%   columns (or rows) of G.  Setting COLIDX to ':' has the same effect as
%   setting it to 1:SIZE(F, 2).
%
%   If F is a row CHEBFUN, then ASSIGNCOLUMNS(F, ROWIDX, G) behaves as
%   described above, except that it assigns the rows of F so that F(ROWIDX,:) =
%   G.
%
%   In both cases, G may also be a numerical vector of the appropriate
%   orientation (column or row).  In this case, G is treated as an array-valued
%   CHEBFUN with a constant value equal to this vector.
%
% See also EXTRACTCOLUMNS, MAT2CELL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% This shouldn't happen:
if ( ~isa(f, 'chebfun') )
    error('CHEBFUN:assigncColumns:notChebfun', 'First input must be a CHEBFUN.')
end

% Number of columns (or rows if f.isTransposed) of f:
numColsF = min(size(f));

% Expand ':' to 1:end:
if ( ~isnumeric(colIdx) && strcmp(colIdx, ':') )
    colIdx = 1:numColsF;
end

% Allow scalar expansion:
if ( isnumeric(g) )
    if ( f.isTransposed )
        g = g.';
        g = chebfun(g, f.domain).';
    else
        g = chebfun(g, f.domain);
    end
end

% Check dimensions of g:
if ( xor(f.isTransposed, g.isTransposed) || (numel(colIdx) ~= min(size(g))) )
    error('CHEBFUN:assignColumns:numCols', ...
        'Subscripted assignment dimension mismatch.')
end

% Trivial case:
if ( (numel(colIdx) == numColsF) && isequal(colIdx, 1:numColsF) )
    % Verify domain:
    if ( ~domainCheck(f, g) )
        error('CHEBFUN:assignColumns:domain', ...
            'Inconsistent domains; domain(f) ~= domain(g).');
    end
    f = g;
    return
end

% Check dimensions of f:
if ( max(colIdx) > numColsF )
    error('CHEBFUN:assignColumns:dims', 'Index exceeds CHEBFUN dimensions.')
end

% Make sure f and g have the same breakpoints:
[f, g] = overlap(f, g);

% Assign the columns of the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = assignColumns(f.funs{k}, colIdx, g.funs{k});
end

% Assign the columns to the impulses:
f.impulses(:,colIdx,:) = g.impulses;

end
