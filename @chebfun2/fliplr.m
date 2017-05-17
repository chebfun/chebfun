function varargout = fliplr(varargin)
%FLIPLR   Flip/reverse a CHEBFUN2 in the x-direction.
%   G = FLIPLR( F ) returns a CHEBFUN2 G with the same domain as F but reversed;
%   that is, G(x,y) = F(a+b-x,y), where the domain is [a, b, c, d].
%
% See also FLIPUD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = fliplr@separableApprox(varargin{:});

end
