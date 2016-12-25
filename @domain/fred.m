function F = fred(k, d, varargin)
%FRED      Fredholm integral operator.
%   This function is deprecated. Use OPERATORBLOCK.FRED.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = linop( operatorBlock.fred(k, double(d), varargin{:}) );

end
