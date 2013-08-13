function f = conj(f)
%CONJ   Complex conjugate of a CHEBFUN.
%   CONJ(F) is the complex conjugate of F.
%
% See also CHEBFUN/REAL, CHEBFUN/IMAG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Conjugate the impulses:
f.impulses = conj(f.impulses);

% Conjugate the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = conj(f.funs{k});
end

end
