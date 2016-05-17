function varargin = conv(varargout)
%CONV   Operation not supported for TRIGTECH objects.
%
%   Use the function CIRCCONV.
%
% See also CIRCCONV.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:TRIGTECH:conv:NotSupported', ...
    'CONV is not supported for trigtech objects. Use the function CIRCCONV.');

end
