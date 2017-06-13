function g = cummin(f)
%CUMMIN   Cumulative minimum of a CHEBFUN.
%   G = CUMMIN(F) is the cumulative minimum of a row or column CHEBFUN F
%   over its domain of definition.
%
% Example:
%
%   f = chebfun('t*sin(t)',[0 60]); plot(f,'b',cummin(f),'r')
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
if ( (rows == Inf && cols > 1) || (cols == Inf && rows > 1) )
    error('CHEBFUN:CHEBFUN:cummin:quasi', ...
        'CUMMIN does not currently support quasimatrices.');
end

% Solve the problem by calling CUMMAX:
g = -cummax(-f);

end