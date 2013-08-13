function f = imag(f)
%IMAG   Complex imaginary part of a CHEBFUN.
%   IMAG(F) is the imaginary part of F.
%
% See also CHEBFUN/REAL.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Take imaginary part of the impulses:
f.impulses = imag(f.impulses);

% Take imaginary part of the funs:
for k = 1:numel(f.funs)
    f.funs{k} = imag(f.funs{k});
end

end
