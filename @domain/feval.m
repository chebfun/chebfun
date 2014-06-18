function F = feval(d, loc, varargin)
%FEVAL     Evaluation functional.
%   This function is deprecated. Use FUNCTIONALBLOCK.FEVAL.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = linop( functionalBlock.feval(loc, d, varargin{:}) );

end
