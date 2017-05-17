function varargout = legcoeffs(varargin) %#ok<STOUT>
%LEGCOEFFS   LEGCOEFFS does not support SINGFUN objects.
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Add support for this.

error('CHEBFUN:SINGFUN:legcoeffs:nosupport', ...
      'LEGCOEFFS does not support SINGFUN objects.')
    
end
