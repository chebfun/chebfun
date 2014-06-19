function varargout = poly(varargin) %#ok<STOUT>
%POLY   POLY does not support SINGFUN objects.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SINGFUN:poly:noSupport', ...
      'POLY does not support SINGFUN objects.')
    
end
