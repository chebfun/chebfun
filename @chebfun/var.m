function out = var(f)
%VAR   Variance of a CHEBFUN.
%   VAR(F) is the variance of the CHEBFUN F.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if  ( ~f(1).isTransposed )
    Y = f - mean(f);
    out = mean(Y.*conj(Y));
else
    out = transpose(var(transpose(f)));
end

end

