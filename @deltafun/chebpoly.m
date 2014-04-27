function out = chebpoly(f, varargin)
%CHEBPOLY   Chebyshev polynomial coefficients of a DELTAFUN.
%   CHEBPOLY(F) returns the Chebyshev coefficients of F.FUNPART.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call CHEBPOLY() of the .FUNPART:
out = chebpoly(f.funPart, varargin{:});

end
