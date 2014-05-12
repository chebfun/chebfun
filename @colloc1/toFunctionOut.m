function f = toFunctionOut(varargin)
%TOFUNCTIONOUT   Convert COLLOC1 discretization to a CHEBFUN. 
%   F = TOFUNCTIONOUT(...) is equivalent to TOFUNCTIONIN(...) for COLLOC1.
%
% See also TOVALUES, TOFUNCTIONIN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = toFunctionIn(varargin{:});

end
