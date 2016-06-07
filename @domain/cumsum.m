function J = cumsum(d, varargin)
%CUMSUM    Indefinite integration operator.
%   This function is deprecated. Use OPERATORBLOCK.CUMSUM.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

J = linop( operatorBlock.cumsum(double(d), varargin{:}) );

end
