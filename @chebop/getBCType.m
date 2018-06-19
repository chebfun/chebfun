function bcType = getBCType(N)
%GETBCTYPE   Detects the type of boundary conditions of a CHEBOP.
%   BCTYPE = GETBCTYPE(N) returns a string identifying the type of boundary
%   conditions used by a CHEBOP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(N.bc) )
    bcType = 'bvp';
elseif ( strcmpi(N.bc, 'periodic') && isempty(N.lbc) && isempty(N.rbc) )
    bcType = 'periodic';
else
    bcType = 'general';
end

end