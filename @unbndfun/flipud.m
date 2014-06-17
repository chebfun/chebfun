function varargout = flipud(varargin) %#ok<STOUT>
%FLIPUD   Flip/reverse an UNBNDFUN object.
%   FLIPUD is not support for FUN objects on unbounded domains.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:UNBNDFUN:flipud:noSupport', ...
    'FLIPUD does not support UNBNDFUN objects.');

end
