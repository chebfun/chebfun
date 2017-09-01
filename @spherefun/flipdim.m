function varargout = flipdim(varargin)
%FLIPDIM   Flip/reverse a SPHEREFUN in a chosen direction.
%   G = FLIPDIM(F, DIM) returns a SPHEREFUN G with the same domain as F but
%   reversed in a direction, i.e., G(x,y)=F(x, c+d-y). If DIM = 2 (default) then
%   G(x,y) = F(x, c+d-y).  Otherwise DIM = 1 and G(x,y) = F(a+b-x, y). The
%   domain of F is [a, b, c, d].
% 
% See also SPHEREFUN/FLIPLR, SPHEREFUN/FLIPUD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = flipdim@separableApprox(varargin{:});

end
