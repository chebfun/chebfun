function V = volt(k, d, varargin)
%VOLT      Volterra integral operator.
%   This function is deprecated. Use OPERATORBLOCK.VOLT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

V = linop( operatorBlock.volt(d, k, varargin{:}) );

end
