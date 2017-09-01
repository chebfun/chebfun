function out = iszero(F, varargin)
%ISZERO   Check if a CHEBFUN is identically zero on its domain.
%   ISZERO(F) returns true if F is identically zero or empty on F.domain and
%   false otherwise. If F is an array-valued CHEBFUN, the a true/false value is
%   returned for each column.

% TODO:  Document the TOL input.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% An empty CHEBFUN is zero:
if ( isempty(F) )
    out = true;
    return
end

out = cell(1, numel(F));
for k = 1:numel(F)
    out{k} = columnIszero(F(k), varargin{:});
end
out = cell2mat(out);
if ( F(1).isTransposed )
    out = out.';
end

end

function out = columnIszero(f, tol)

% Choose a tolerance:
if ( nargin < 2 )
    tol = vscale(f)*eps;
end

% pointValues:
out = all(f.pointValues <= tol);

% Loop over each of the FUNs:
% TODO:  We don't use TOL here?
k = 0;
while ( k < numel(f.funs) && any(out) )
    k = k + 1;
    out = out & iszero(f.funs{k});
end

end
