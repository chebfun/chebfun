function f = addBreaksForCumSum(f)
%ADDBREAKSFORCUMSUM   Introduce a new breakpoint to facilitate the computation 
%   of the indefinite integral of a BNDFUN.
%
%   G = ADDBREAKSFORCUMSUM(F) returns a cell of two BNDFUN when F is made of 
%   SINGFUN whose exponents are both non-trivial. Note that this function is 
%   expected to be called by CUMSUM only.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Is F made of SINGFUN?
if ( ~issing(f) )
    return
end

% Non-trivial exponents at both endpoints?
if ( any(~get(f, 'exponents') ) )
    return
end

% Breakpoints are the mid-point of each domain which needs to be broken:
breakPoint = f.mapping.for(0);

% Introduce new break points using RESTRICT.
f = restrict(f, [f.domain(1) breakPoint f.domain(2)]);

end