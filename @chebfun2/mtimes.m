function varargout = mtimes(varargin)
%*	   Pointwise multiplication for CHEBFUN2 objects.
%
%   c*F or F*c multiplies a CHEBFUN2 F by a scalar c.
%
%   F*G computes the integral of F(s,y)G(x,s) over s, and this is the continuous
%   analogue of matrix-matrix multiplication.
%
% See also TIMES.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mtimes@separableApprox(varargin{:});

end
