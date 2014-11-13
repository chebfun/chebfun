function f = toFunctionIn(disc, values)
%TOFUNCTIONIN   Convert TRIGCOLLOC discretization to a CHEBFUN. 
%   TOFUNCTIONIN(DISC, VALUES) converts the values of a TRIGCOLLOC-discretized
%   function to a CHEBFUN. 
%
%   If VALUES is matrix valued, the output is an array-valued CHEBFUN, where
%   each column of the CHEBFUN corresponds to a column of the input.
%
% See also TOVALUES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Convert the VALUES matrix into a CHEBFUN on the appropriate domain
% (potentially an array-valued CHEBFUN).
f = chebfun(values, disc.domain, 'periodic');

end