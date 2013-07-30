function out = var(f)
%VAR   Variance of a chebfun.
%   VAR(F) is the variance of the chebfun F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ~( f.isTransposed )
    Y = f - mean(f);
    out = mean(Y.*conj(Y));
else
    out = transpose(var(transpose(f)));
end

end

