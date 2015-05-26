function F = diff(d, varargin)
%DIFF      Differentiation operator.
%   This function is deprecated. Use OPERATORBLOCK.DIFF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = linop( operatorBlock.diff(double(d), varargin{:}) );

end
