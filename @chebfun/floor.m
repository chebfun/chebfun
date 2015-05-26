function f = floor(f)
%FLOOR   Pointwise floor function of a CHEBFUN.
%   G = FLOOR(F) returns the CHEBFUN G such that G(X) = FLOOR(F(x)) for each x
%   in F.domain. 
%
%   If F is complex, then the G = FLOOR(REAL(F)) + 1i*FLOOR(IMAG(F)).
%
% See also CEIL, ROUND, FIX.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with the trivial empty case:
if ( isempty(f) )
    return
end

% Deal with unbounded functions:
if ( ~isfinite(f) )
    error('CHEBFUN:CHEBFUN:floor:inf', ...
        'Floor is not defined for functions which diverge to infinity.');
end

% Deal with complex-valued functions:
if ( ~isreal(f) )
    if ( isreal(1i*f) )
        f = 1i*floor(imag(f));
    else
        f = floor(real(f)) + 1i*floor(imag(f));
    end
    return
end

% Find all the integer crossings for f:
mmvals = minandmax(f); % [TODO]: Only need a good bound?
minf = min(mmvals(1, :));
maxf = max(mmvals(2, :));
range = floor([minf, maxf]);
for k = (range(1)+1):range(2)
    f = addBreaksAtRoots(f - k) + k;
end

% Loop over the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = floor(f.funs{k});
end

% Floor the pointValues:
f.pointValues = floor(f.pointValues);

% Simplify the result:
f = merge(f);

end
