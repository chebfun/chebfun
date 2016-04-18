function varargout = mldivide(varargin)
%\   Left divide for SPHEREFUN2 objects.
%
%    Left divide for a SPHEREFUN2. Only allowed to divide by scalars.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = mldivide@separableApprox(varargin{:});
end
