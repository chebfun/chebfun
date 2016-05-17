function f = fix(f)
%FIX   Round a CHEBFUN pointwise toward zero.
%   G = FIX(F) returns the CHEBFUN G such that G(x) = FIX(F(x)) for each x in
%   F.domain.
%
%   If F is complex, then the G = FIX(REAL(F)) + 1i*FIX(IMAG(F)).
%
% See also ROUND, CEIL, FLOOR.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with the trivial empty case:
if ( isempty(f) )
    return
end

% Deal with unbounded functions:
if ( ~isfinite(f) )
    error('CHEBFUN:CHEBFUN:fix:inf', ...
        'fix() is not defined for functions which diverge to infinity.');
end

% Deal with complex-valued functions:
if ( ~isreal(f) )
    if ( isreal(1i*f) )
        f = 1i*fix(imag(f));
    else
        f = fix(real(f)) + 1i*fix(imag(f));
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
    f.funs{k} = fix(f.funs{k});
end

% Fix the pointValues:
f.pointValues = fix(f.pointValues);

% Simplify the result:
f = merge(f);

end
