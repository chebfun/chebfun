function f = round(f)
%ROUND   Round a CHEBFUN pointwise to nearest integer.
%   G = ROUND(F) returns the CHEBFUN G such that G(x) = ROUND(F(x)) for each x
%   in F.domain.
%
%   If F is complex, then the G = ROUND(REAL(F)) + 1i*ROUND(IMAG(F)).
%
% See also FIX, FLOOR, CEIL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with the trivial empty case:
if ( isempty(f) )
    return
end

% Deal with unbounded functions:
if ( ~isfinite(f) )
    error('CHEBFUN:CHEBFUN:round:inf', ...
        'round() is not defined for functions which diverge to infinity.');
end

% Deal with complex-valued functions:
if ( ~isreal(f) )
    if ( isreal(1i*f) )
        f = 1i*round(imag(f));
    else
        f = round(real(f)) + 1i*round(imag(f));
    end
    return
end

% Find all the integer-plus-0.5 crossings for f:
mmvals = minandmax(f); % [TODO]: Only need a good bound?
minf = min(mmvals(1, :));
maxf = max(mmvals(2, :));
range = floor([minf, maxf]);
for k = (range(1) + 1):(range(2) + 1)
    f = addBreaksAtRoots(f - k + .5) + k - .5;
end

% Loop over the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = round(f.funs{k});
end

% Round the pointValues:
f.pointValues = round(f.pointValues);

% Simplify the result:
f = merge(f);

end

