function g = nextpow2(f, pref)
%NEXTPOW2   Base 2 power of a CHEBFUN.
%   P = NEXTPOW2(N) returns the first P such that 2.^P >= abs(N). It is often
%   useful for finding the nearest power of two sequence length for FFT
%   operations.
%
% See also LOG2, POW2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information. 

% Trivial case:
if ( isempty(f) )
    g = f;
    return
elseif ( ~isreal(f) )
    error('CHEBFUN:nextpow2:complex', 'Input to NEXTPOW() must be real.'); 
end

% Grab preferences:
if ( nargin < 2 )
    pref = chebpref();
end

% Compute the absolute value of f:
absf = abs(f);

% Result will have breaks where |f| is a power of 2.
mm = minandmax(absf);
pows = 2.^(nextpow2(mm(1,:)):nextpow2(mm(2,:)));
r = f.domain.';
for k = 1:numel(pows)
    rk = roots(absf - pows(k));
    r = [r ; rk(:)];
end
r = unique(r);
r(isnan(r)) = [];

% Restrict f to these points:
f = restrict(f, r.');

% We need to extrapolate:
pref.extrapolate = true;

% Call COMPOSE():
g = compose(f, @(n) nextpow2(n), [], pref);

% Attempt to remove unnecessary breaks:
g = merge(g); % TODO: Could we work this out in advance?

end
