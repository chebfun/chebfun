function varargout = flipud(varargin)
%FLIPUD   Flip/reverse a CHEBFUN2 in the y-direction.
%   G = FLIPUD(F) returns a CHEBFUN2 G with the same domain as F but reversed;
%   that is, G(x,y) = F(x, c+d-y), where the domain is [a, b, c, d].
%
% See also FLIPLR, FLIPDIM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = flipud@separableApprox(varargin{:});

end
