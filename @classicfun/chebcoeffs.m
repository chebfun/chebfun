function out = chebcoeffs(f, varargin)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a CLASSICFUN.
%   CHEBCOEFFS(F) returns the Chebyshev coefficients of F.ONEFUN.
%
% See also LEGPOLY FOURCOEFFS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call CHEBCOEFFS() of the .ONEFUN:
out = chebcoeffs(f.onefun, varargin{:});

end
