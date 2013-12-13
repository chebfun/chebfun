function isf = isfinite(f)
%ISFINITE   Test if a CHEBFUN is bounded.
%   ISFINITE(F) returns FALSE if F has any infinite values and TRUE otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for infinite breakpoint values:
if ( any(any(~isfinite(f.pointValues))) )
    isf = false;
    return
end

% Check for any infinite FUNs.
for k = 1:numel(f.funs)
    if ( ~isfinite(f.funs{k}) )
        isf = false;
        return
    end
end

isf = true;

end
