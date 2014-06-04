function spy(N, varargin)
%SPY    Visualize a linear CHEBOP.
%   SPY(A) creates a picture of the nonzero pattern of the default
%   discretization of the linear CHEBOP A. See LINOP/SPY for details.
%
% See also LINOP/SPY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Try to create a LINOP. An error will be thrown if n is not linear.
L = linop(N);

% SPY plot of the LINOP:
spy(L, varargin{:});

end