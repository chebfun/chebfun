function F = diff(d, varargin)
%DIFF      Differentiation operator.
%   This function is deprecated. Use OPERATORBLOCK.DIFF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

F = linop( operatorBlock.diff(d, varargin{:}) );

end
