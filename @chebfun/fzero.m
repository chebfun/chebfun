function varargin = fzero(varargout) %#ok<STOUT,INUSD>
%FZERO   Deprecated function.
%   FZERO(F, ...) is no longer supported. Use ROOTS(F).
%
% See also ROOTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:CHEBFUN:roots:deprecated', ...
    'FZERO(F) for CHEBFUN input F is no longer supported. Use ROOTS(F).');

end