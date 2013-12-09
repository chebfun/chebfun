function out = isnan(f)
%ISNAN   Test if a CHEBFUN is NaN.
%   ISNAN(F) returns TRUE if F has any NaN values and FALSE otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = false;

% Empty CHEBFUNs are not NaN.
if ( isempty(f) )
    return;
end

% Check for NaN impulse values.
if ( any(isnan(f.impulses(:))) )
    out = true;
    return;
end

% Check for NaN FUNs.
numCols = size(f.funs{1}, 2);
for k = 1:numCols
    if ( isnan(f.funs{k}) )
        out = true;
        return;
    end
end

end
