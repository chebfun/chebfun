function f = imag(f)
%IMAG   Complex imaginary part.
%   IMAG(F) is the imaginary part of F.
%
% See also REAL.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Take imaginary part of the impulses:
f.impulses = imag(f.impulses);

% Take imaginary part of the funs:
for k = 1:numel(f.funs)
    f.funs{k} = imag(f.funs{k});
end

end