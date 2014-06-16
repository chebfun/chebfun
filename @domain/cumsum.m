function J = cumsum(d, varargin)
%CUMSUM    Indefinite integration operator.
%   This function is deprecated. Use OPERATORBLOCK.CUMSUM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

J = linop( operatorBlock.cumsum(d, varargin{:}) );

end
