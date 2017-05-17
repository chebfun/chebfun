function varargout = minus(varargin)
%-   Subtraction of two SPHEREFUN objects.
% 
%   F - G subtracts G from F, where F and G are SPHEREFUN objects or scalars.
%
% See also SPHEREFUN/PLUS, SPHEREFUN/UMINUS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = minus@separableApprox(varargin{:});

end
