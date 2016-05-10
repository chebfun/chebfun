function varargout = mldivide(varargin)
%\      Left divide for CHEBFUN2 objects.
%
% Left divide for a CHEBFUN2. Only allowed to divide by scalars.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mldivide@separableApprox(varargin{:});

end
