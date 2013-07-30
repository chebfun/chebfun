function varargout = flipud(varargin) %#ok<STOUT>
%FLIPUD   Flip/reverse an UNBNDFUN object.
%   FLIPUD is not support for FUN objects on unbounded domains.

error('CHEBFUN:UNBNDFUN:flipud:nosupport', ...
    'FLIPUD does not support UNBNDFUN objects.');

end