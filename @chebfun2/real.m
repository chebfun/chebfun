function varargout = real(varargin)
%REAL      Real part of a CHEBFUN2.
%
% See also IMAG.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = real@separableApprox(varargin{:});

end
