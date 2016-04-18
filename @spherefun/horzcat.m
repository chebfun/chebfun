function varargout = horzcat(varargin)
%HORZCAT Horizontal concatenation of SPHEREFUN2 objects.
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a SPHEREFUN2 user.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = horzcat@separableApprox(varargin{:});
end
