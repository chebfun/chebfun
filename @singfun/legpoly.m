function varargout = legpoly(varargin) %#ok<STOUT>
%LEGPOLY   LEGPOLY does not support SINGFUN objects.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SINGFUN:legpoly:nosupport', ...
      'LEGPOLY does not support SINGFUN objects.')
    
end