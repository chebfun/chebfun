function varargin = quad(varargout) %#ok<STOUT,INUSD>
%QUAD   Deprecated function.
%   I = QUAD(F, ...) is no longer supported. Use I = INTEGRAL(F, ...).
%
% See also INTEGRAL, SUM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:CHEBFUN:quad:deprecated', ...
    'QUAD(F) for CHEBFUN input F is no longer supported. Use INTEGRAL(F).');

end