function f = assignColumns(f, colIdx, g)
%ASSIGNCOLUMNS   Assign columns (or rows) of an array-valued CHEBFUN.
%   F = ASSIGNCOLUMNS(F, COLIDX, G), if F is a column CHEBFUN, assigns the
%   columns of the CHEBFUN G (or rows, if G is a row CHEBFUN) to the columns
%   specified by the vector COLIDX so that F(:,COLIDX) = G. COLIDX need not be
%   increasing in order or unique but must contain only integers in the range
%   [1, M] (where F has M columns) and have a length equal to the number of
%   columns (or rows) of G. Setting COLIDX to ':' has the same effect as setting
%   it to 1:SIZE(F, 2). If G is empty, then columns of F are removed.
%
%   If F is a row CHEBFUN, then ASSIGNCOLUMNS(F, ROWIDX, G) behaves as described
%   above, except that it assigns the rows of F so that F(ROWIDX,:) = G.
%
%   In both cases, G may also be a numerical vector of the appropriate
%   orientation (column or row).  In this case, G is treated as an array-valued
%   CHEBFUN with a constant value equal to this vector.
%
% See also EXTRACTCOLUMNS, MAT2CELL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% This shouldn't happen:
if ( ~isa(f, 'chebfun') )
    error('CHEBFUN:CHEBFUN:assignColumns:notChebfun', ...
        'First input must be a CHEBFUN.')
end

% Number of columns (or rows if f.isTransposed) of f:
numColsF = numColumns(f);

% Expand ':' to 1:end:
if ( ~isnumeric(colIdx) && strcmp(colIdx, ':') )
    colIdx = 1:numColsF;
end

% g is empty - Remove columns!
if ( isempty(g) )
    if ( numel(f) > 1 )
        f(colIdx) = [];
    else
        for k = 1:numel(f.funs)
            f.funs{k} = assignColumns(f.funs{k}, colIdx, []);
        end
        f.pointValues(:,colIdx) = [];
    end
    if ( numel(f) == 0 )
        % Create an empty CHEBFUN if we have removed all columns:
        f = chebfun();
    end
    return
end

% Allow scalar expansion:
if ( isnumeric(g) )
    if ( f(1).isTransposed )
        g = g.';
        g = chebfun(g, domain(f, 'ends')).';
    else
        g = chebfun(g, domain(f, 'ends'));
    end
end

% Check dimensions of g:
if ( xor(f(1).isTransposed, g(1).isTransposed) || (numel(colIdx) ~= numColumns(g)) )
    error('CHEBFUN:CHEBFUN:assignColumns:numCols', ...
        'Subscripted assignment dimension mismatch.')
end

% Trivial case:
if ( (numel(colIdx) == numColsF) && isequal(colIdx, 1:numColsF) )
    % Verify domain:
    if ( ~domainCheck(f, g) )
        error('CHEBFUN:CHEBFUN:assignColumns:domain', ...
            'Inconsistent domains; domain(f) ~= domain(g).');
    end
    f = g;
    return
end

% Check dimensions of f:
if ( max(colIdx) > numColsF )
%     error('CHEBFUN:CHEBFUN:assignColumns:dims', 'Index exceeds CHEBFUN dimensions.')

    if ( isempty(f) )
        dom = g.domain;
    else
        dom = f.domain;
    end
    
    % Pad with zeros so that f has sufficiently many columns:
    z = chebfun(zeros(1, max(colIdx) - numColsF), dom);
    f = [f z];

end

if ( numel(f) == 1 && numel(g) == 1)
    % Array-valued CHEBFUN case:

    % Make sure f and g have the same breakpoints:
    [f, g] = overlap(f, g);
    % Assign the columns of the FUNs:
    for k = 1:numel(f.funs)
        f.funs{k} = assignColumns(f.funs{k}, colIdx, g.funs{k});
    end
    % Assign the columns to the impulses:
    f.pointValues(:, colIdx, :) = g.pointValues;
else
    % Quasimatrix case:
    
    % Ensure f is a quasimatrix:
    f = cheb2quasi(f);
    g = cheb2cell(g);
    for k = 1:numel(colIdx)
        f(colIdx(k)) = g{k};
    end
end
