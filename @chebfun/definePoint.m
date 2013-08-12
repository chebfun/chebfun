function f = definePoint(f, s, v)
% DEFINE Supply a new definition for a chebfun at a point / set of points.
%
%   F = DEFINEPOINT(F, S, V) redefines the CHEBFUN F to take the values V at the
%   points S in F.domain. If F is a scalar-valued CHEBFUN, then S and V should
%   vectors of the same length. If F is an array-valued CHEBFUN, then V should
%   be a matrix of size length(S)*size(F, 2).
%
%   An equivalent syntax is F(S) = V.
%
% See also CHEBFUN/SUBSASGN, CHEBFUN/RESTRICT, CHEBFUN/DEFINEINTERVAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%[TODO]: Documentation (particularly to describe dimensions and array-valued
%        Chebfun objects.

% Trival empty case:
if ( isempty(s) )
    return
end

% Some error checking:
if ( isempty(v) )
    error('CHEBFUN:subsasgn:conversion',...
            'Cannot assign empty values to points.')
elseif ( ~isa(v, 'numeric') )
    error('CHEBFUN:subsasgn:conversion',...
            ['Conversion to numeric from ', class(v), ' is not possible.'])
elseif ( min(s) < f.domain(1) || max(s) > f.domain(end) )
    error('CHEBFUN:subsasgn:outbounds',...
        'Cannot introduce points outside the domain.')
end

% Deal with row chebfuns:
dim = 2 - f.isTransposed;
% Find the number of columns (or rows in a row chebfun):
numCols = size(f, dim);

% Expand v if required:
if ( length(v) == 1 )
    v = repmat(v, numel(s), numCols);
elseif ( length(numCols) == 1 && min(size(v)) == 1 && length(v) == length(s) )
    v = v(:);
elseif ( numel(s) ~= size(v, 1) || numCols ~= size(v, 2) )
    error('CHEBFUN:subsasgn:dimensions',...
            'Subscripted assignment dimension mismatch.')
end

% Restrict f to the new subdomains:
ss = [f.domain(1) ; s(:) ; f.domain(end)]';
f = restrict(f, ss);

% Assign the values in v to the new impulses;
[mem, loc] = ismember(s, f.domain);
f.impulses(loc,:,1) = v(mem, :);

end