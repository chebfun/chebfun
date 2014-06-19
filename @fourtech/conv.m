function varargin = conv(varargout)
%CONV   Operation not supported for FOURTECH objects.
%
%   Use the function CIRCCONV.
%
% See also CIRCCONV.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:FOURTECH:conv:NotSupported', ...
    'CONV is not supported for fourtech objects. Use the function CIRCCONV.');

end