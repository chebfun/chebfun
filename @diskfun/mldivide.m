function varargout = mldivide(varargin)
%\   Left divide for DISKFUN objects.
%   MLDIVIDE(F, G) computes left dividison for DISKFUN objects F and G.
%   Only allowed to divide by scalars.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mldivide@separableApprox(varargin{:});

end