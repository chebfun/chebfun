function g = cummin(f)
%CUMMAX   Cumulative minimum of a CHEBFUN.
%   G = CUMMIN(F) 
%
% See also CUMMAX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case:
if ( isempty(f) )
    return
end

% Quasimatrices are not supported:
[rows, cols] = size(f);
if ( rows == Inf || cols == Inf )
    error('CHEBFUN:CHEBFUN:cummin:quasi', ...
        'CUMMIN does not currently support quasimatrices.');
end

% Call CUMMAX:
g = -cummax(-f);

end