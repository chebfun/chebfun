function varargout = mldivide(varargin)
%\   Left divide for DISKFUN objects.
%
%    Left divide for a DISKFUN. Only allowed to divide by scalars.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mldivide@separableApprox(varargin{:});

end
