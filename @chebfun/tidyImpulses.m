function f = tidyImpulses(f)
%TIDYIMPULSES   Remove all-zero layers of higher-order impulses.
%   G = TIDYIMPULSES(F) returns a CHEBFUN G which is identical to the CHEBFUN
%   F except that any trailing all-zero layers of higher-order impulses are
%   stripped away.

% [TODO]: This function should no longer be required since the DELTAFUN class
% handles all delta functions.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

for n = size(f.pointValues, 3):-1:2
    if ( all(~any(f.pointValues(:, :, n))) )
        f.pointValues(:,:,n) = [];
    else
        break
    end
end

end
