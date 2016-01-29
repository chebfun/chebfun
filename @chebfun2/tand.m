function varargout = tand(varargin)
%TAND  Tangent of a CHEBFUN2 (in degrees)

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = tand@separableApprox(varargin{:});

end
