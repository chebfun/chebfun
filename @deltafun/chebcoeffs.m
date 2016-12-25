function out = chebcoeffs(f, varargin)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a DELTAFUN.
%   CHEBCOEFFS(F) returns the Chebyshev coefficients of F.FUNPART.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call CHEBCOEFFS() of the .FUNPART:
out = chebcoeffs(f.funPart, varargin{:});

end
