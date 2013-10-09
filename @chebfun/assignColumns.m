function f = assignColumns(f, colIdx, g)
%ASSIGNCOLUMNS   Assign columns (or rows) of an array-valued CHEBFUN.
%   F = ASSIGNCOLUMNS(F, COLIDX, G), if F is a column CHEBFUN, assigns the
%   columns of G (or rows, if G is a row CHEBFUN) to the columns specified by
%   the vector COLIDX so that F(:,COLIDX) = G. COLIDX need not be increasing in
%   order or unique but must contain only integers in the range [1, M] (where F
%   has M columns) and have a length equal to the number of columns (or rows)
%   of G.
%
%   If F is a row CHEBFUN, then ASSIGNCOLUMNS(F, ROWIDX) behaves as described
%   above, except that it assigns the rows of F so that F(ROWIDX,:) = G.
%
% See also EXTRACTCOLUMNS, MAT2CELL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Number of columns (or rows if f.isTransposed) of f:
numColsF = min(size(f));
numColsG = min(size(g));

% Trivial cases:
if ( numel(colIdx) ~= min(size(g)) )
    error('CHEBFUN:ASSIGNCOLUMNS:numCols', 'Index exceeds CHEBFUN dimensions.')
end
if ( ~isnumeric(colIdx) && strcmp(colIdx, ':') )
    f = g;
elseif ( (numel(colIdx) == numel(numColsF)) && all(colIdx == 1:numColsF) )
    f = g;
elseif ( max(colIdx) > numColsF )
    error('CHEBFUN:subsref:dimensions', 'Index exceeds CHEBFUN dimensions.')    
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
