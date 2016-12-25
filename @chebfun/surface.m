function varargout = surface(varargin)
%SURFACE   Surface plot for CHEBFUN objects.
%   SURFACE(...) is the same as SURF(...).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = surf(varargin{:});

end
