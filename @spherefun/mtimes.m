function varargout = mtimes(varargin)
%*	   Pointwise multiplication for SPHEREFUN objects.
%
%   c*F or F*c multiplies a SPHEREFUN F by a scalar c.
%
%   F*G computes the integral of F(s,t)G(l,s) over s, and this is the 
%   continuous analogue of matrix-matrix multiplication.
%
% See also SPHEREFUN/TIMES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mtimes@separableApprox(varargin{:});

end
