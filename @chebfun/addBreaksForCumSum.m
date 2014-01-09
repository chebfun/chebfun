function f = addBreaksForCumSum(f)
%NEEDTOBREAK   Introduce new breakpoints to facilitate the computation of the 
%   indefinite integral of a CHEBFUN.
%
%   G = ADDBREAKSFORCUMSUM(F) returns a CHEBFUN which may have new breakpoint which is 
%   necessary when F has at least one FUN which is made of SINGFUN whose 
%   exponents are both non-trivial. Note that this function is expected to be 
%   called in CUMSUM only. Users are not encouraged to use this function.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Index for FUNs which are made of SINGFUN:
isSing = cellfun(@(f) isa(f.onefun, 'singfun'), f.funs);

% Index for FUNs which are SINGFUN with non-trivial exponents at both endpoints:
nonZeroExps = cellfun( @(f) any(f.onefun.exponents), f.funs(isSing) );

% Index of FUNs which need to be broken:
ind = isSing;
ind(isSing) = nonZeroExps;

% Breakpoints are the mid-point of each domain which needs to be broken:
dom = f.domain;
breakPoints = (dom(1:end-1) + dom(2:end))/2;
breakPoints = breakPoints(ind);

% Introduce new break points using RESTRICT.
if ( ~isempty(breakPoints) )
    domNew = sort([f.domain, breakPoints]);
    f = restrict(f, domNew);
end

end