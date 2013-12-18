function f = sqrt(f)
%SQRT   Square root of a CHEBFUN.
%   SQRT(F) returns the square root of a CHEBFUN F.
%
% See also POWER.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case: (f is empty)
if ( isempty(f) )
    return
end

% Add breaks at the appropriate roots of f:
if ( isreal(f) )
    f = addBreaksAtRoots(f);
else
    f = addBreaksAtRoots(f, 'imag');
end

% Loop over each FUN and call SQRT@BNDFUN on each of the FUNs:
numFuns = numel(f.funs);
for k = 1:numFuns
    f.funs{k} = sqrt(f.funs{k});
end

% Take the absolute value of the impulses in the first row:
f.impulses = abs(f.impulses(:,:,1));

% [TODO]: Is simplify needed?
f = simplify(f);

end