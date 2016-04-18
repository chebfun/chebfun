function varargout = minus(varargin)
%-   Subtraction of two SPHEREFUN2 objects.
% 
%   F - G subtracts G from F, where F and G are SPHEREFUN2 objects or scalars.
%
% See also PLUS, UMINUS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = minus@separableApprox(varargin{:});
end
