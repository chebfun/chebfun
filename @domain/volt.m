function V = volt(k, d, varargin)
%VOLT      Volterra integral operator.
%   This function is deprecated. Use OPERATORBLOCK.VOLT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

V = linop( operatorBlock.volt(d, k, varargin{:}) );

end