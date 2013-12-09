function f = definePoint(f, s, v)
%DEFINEPOINT   Supply new definition for a CHEBFUN at a point or set of points.
%   F = DEFINEPOINT(F, S, V) redefines the CHEBFUN F to take the values V at the
%   points S in F.DOMAIN. If F is a scalar-valued CHEBFUN, then S and V should
%   vectors of the same length. If F is an array-valued CHEBFUN, then V should
%   be a matrix of size LENGTH(S)*SIZE(F, 2).
%
%   An equivalent syntax is F(S) = V.
%
% See also SUBSASGN, RESTRICT, DEFINEINTERVAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Trival empty case:
if ( isempty(s) )
    return
end

% Some error checking:
if ( isempty(v) )
    error('CHEBFUN:definePoint:conversion',...
            'Cannot assign empty values to points.')
elseif ( ~isa(v, 'numeric') )
    error('CHEBFUN:definePoint:conversion',...
            ['Conversion to numeric from ', class(v), ' is not possible.'])
elseif ( (min(s) < f.domain(1)) || (max(s) > f.domain(end)) )
    error('CHEBFUN:definePoint:outbounds',...
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
    error('CHEBFUN:definePoint:dimensions',...
            'Subscripted assignment dimension mismatch.')
end

% Restrict f to the new subdomains:
snew = [f.domain(1) ; s(:) ; f.domain(end)]';
f = restrict(f, snew);

% Assign the values in v to the new impulses;
[mem, loc] = ismember(s, f.domain);
f.impulses(loc,:,1) = v(mem, :);

end
