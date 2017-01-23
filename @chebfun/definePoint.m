function f = definePoint(f, s, v)
%DEFINEPOINT   Supply new definition for a CHEBFUN at a point or set of points.
%   F = DEFINEPOINT(F, S, V) redefines the CHEBFUN F to take the values V at the
%   points S in F.DOMAIN. If F is a scalar-valued CHEBFUN, then S and V should
%   be vectors of the same length. If F is an array-valued CHEBFUN or a
%   quasimatrix, then V should be a matrix of size LENGTH(S)*NUMCOLUMNS(F).
%
%   An equivalent syntax is F(S) = V.
%
% See also SUBSASGN, RESTRICT, DEFINEINTERVAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trival empty case:
if ( isempty(s) )
    return
end


if ( numel(f) == 1 )
    % Array-valued CHEBFUN case:
    
    f = columnDefinePoint(f, s, v);
    
else
    % Quasimatrix case:
    
    % Expand v if required:
    if ( size(v, 2 - f(1).isTransposed) == 1 )
        v = repmat(v, 1, numColumns(f));
    end
    if ( f(1).isTransposed )
        v = v.';
    end
    % Loop over the columns:
    for k = numel(f)
        f(k) = columnDefinePoint(f(k), s, v(:,k));
    end
end

end

function f = columnDefinePoint(f, s, v)

% Some error checking:
if ( isempty(v) )
    error('CHEBFUN:CHEBFUN:definePoint:columnDefinePoint:conversion',...
            'Cannot assign empty values to points.')
elseif ( ~isa(v, 'numeric') )
    error('CHEBFUN:CHEBFUN:definePoint:columnDefinePoint:conversion',...
            ['Conversion to numeric from ', class(v), ' is not possible.'])
elseif ( (min(s) < f.domain(1)) || (max(s) > f.domain(end)) )
    error('CHEBFUN:CHEBFUN:definePoint:columnDefinePoint:outbounds',...
        'Cannot introduce points outside the domain.')
end

% Deal with row CHEBFUN objects:
if ( ~f.isTransposed )
    dim = 2;
else
    dim = 1;
end

% Find the number of columns (or rows in a row chebfun):
numCols = size(f, dim);

% Expand v if required:
if ( length(v) == 1 )
    v = repmat(v, numel(s), numCols);
elseif ( (size(v, 2) == numCols) && (size(v, 1) == 1) )
    v = repmat(v, numel(s), 1);
elseif ( (numCols == 1) && (min(size(v)) == 1) && (length(v) == length(s)) )
    v = v(:);
elseif ( (numel(s) ~= size(v, 1)) || (numCols ~= size(v, 2)) )
    error('CHEBFUN:CHEBFUN:definePoint:columnDefinePoint:dimensions',...
            'Subscripted assignment dimension mismatch.')
end

% Restrict f to the new subdomains:
sNew = [f.domain(1) ; s(:) ; f.domain(end)]';
sNew = unique(sNew);
f = restrict(f, sNew);

% Assign the values in v to the new pointValues;
[mem, loc] = ismember(s, f.domain);
f.pointValues(loc,:) = v(mem, :);

end
