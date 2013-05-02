function varargout = poly(varargin) %#ok<STOUT>
%POLY   Polynomial coefficients of an UNBNDFUN object.
%   POLY(F) is not supported for FUN objects on unbounded domains.

error('CHEBFUN:unbndfun:poly:nosupport', ...
    'POLY does not support UNBNDFUN objects.');

end