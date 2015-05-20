function L = addbc(L, varargin)
%ADDBC   Synonym for ADDCONSTRAINT.
%
% See also LINOP.ADDCONSTRAINT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

L = L.addConstraint(varargin{:});

end
