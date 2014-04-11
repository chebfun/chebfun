function varargout = poly(varargin) %#ok<STOUT>
%POLY   Polynomial coefficients of an UNBNDFUN object.
%   POLY(F) is not supported for FUN objects on unbounded domains.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

error('CHEBFUN:unbndfun:poly:nosupport', ...
    'POLY does not support UNBNDFUN objects.');

end
