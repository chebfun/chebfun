function varargout = real(varargin)
%REAL      Real part of a SPHEREFUN.
%
% See also SPHEREFUN/IMAG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = real@separableApprox(varargin{:});
end
