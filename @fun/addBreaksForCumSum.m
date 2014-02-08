function f = addBreaksForCumSum(f)
%ADDBREAKSFORCUMSUM   Introduce a new breakpoint to facilitate the computation 
%   of the indefinite integral of a FUN.
%
%   G = ADDBREAKSFORCUMSUM(F) returns a cell of two FUNs when F is made of 
%   SINGFUN whose exponents are both non-trivial. Note that this function is 
%   expected to be called by CUMSUM only.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% If F is not a SINGFUN or a SINGFUN with only one nonzero exponent, then return
% without doing anything:

if ( ~issing(f) || any(~get(f, 'exponents')) )
    return
end

% Breakpoints are the mid-point of each domain which needs to be broken:
breakPoint = f.mapping.for(0);

% Introduce new break points using RESTRICT.
f = restrict(f, [f.domain(1) breakPoint f.domain(2)]);

end