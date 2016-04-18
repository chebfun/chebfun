function varargout = tan(varargin)
%TAN   Tangent of a SPHEREFUN2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = tan@separableApprox(varargin{:});
end
