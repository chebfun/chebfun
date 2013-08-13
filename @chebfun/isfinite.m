function isf = isfinite(f)
%ISFINITE   Test if a CHEBFUN is bounded.
%   ISFINITE(F) returns FALSE if F has any infinite values and TRUE otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for infinite breakpoint values or higher-order impulses.
if ( any(any(~isfinite(f.impulses(:,:,1)))) || (size(f.impulses, 3) > 1) )
    isf = false;
    return
end

% Check for any infinite funs.
for k = 1:numel(f.funs)
    if ( ~isfinite(f.funs{k}) )
        isf = false;
        return
    end
end

isf = true;

end
