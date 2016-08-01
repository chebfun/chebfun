function varargout = fliplr(varargin)
%FLIPLR   Flip/reverse a DISKFUN in the longitude-direction.
%   G = FLIPLR( F ) returns a DISKFUN G with the same domain as F but reversed;
%   that is, G(x,y) = F(a+b-x,y), where the domain is [a, b, c, d].
%
% See also DISKFUN/FLIPUD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = fliplr@separableApprox(varargin{:});

end
