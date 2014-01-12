function g = inv2(f, varargin)
%INV2   Invert a CHEBFUN.
%   FINV = INV2(F, ...) is a shortcut for FINV = INV(F, ..., 'ALGORITHM',
%   'NEWTON').  For other valid input parameters, se INV.
%
%   Example:
%      f = chebfun(@(x) tanh(7*x)./tanh(7) + 1, [-.5, .5]);
%      g = inv2(f, 'splitting', 'off', 'rangecheck', 'off', 'monocheck', 'off');
%
%   NB:  This function is experimental and slow!  Use of INV with the default
%   'ROOTS' algorithm may be the better choice for piecewise functions, whereas
%   the 'NEWTON' algorithm used by INV2 is good for smooth functions.
%
% See also INV.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = inv(f, varargin{:}, 'algorithm', 'newton');

end
