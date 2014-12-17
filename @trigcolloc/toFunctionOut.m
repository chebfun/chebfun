function f = toFunctionOut(varargin)
%TOFUNCTIONOUT   Convert TRIGCOLLOC discretization to a CHEBFUN. 
%   F = TOFUNCTIONOUT(...) is equivalent to TOFUNCTIONIN(...) for TRIGCOLLOC.
%
% See also TOVALUES, TOFUNCTIONIN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = toFunctionIn(varargin{:});

end