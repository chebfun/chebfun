function varargout = restrict(varargin)
% RESTRICT  Restrict the domain of a CHEBFUN2.
%
% F = RESTRICT(F, DOM) returns a CHEBFUN2 on the domain DOM that approximates F
% F on that domain.  DOM should be a vector of length 4 giving the coordinates
% of the corners.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = restrict@separableApprox(varargin{:});

end
