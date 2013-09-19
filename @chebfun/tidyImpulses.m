function f = tidyImpulses(f)
%TIDYIMPULSES   Remove all-zero layers of higher-order impulses.
%   G = TIDYIMPULSES(F) returns a CHEBFUN G which is identical to the CHEBFUN
%   F except that any trailing all-zero layers of higher-order impulses are
%   stripped away.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

for n = size(f.impulses, 3):-1:1
    if ( all(~any(f.impulses(:, :, n))) )
        f.impulses(:,:,n) = [];
    else
        break
    end
end

end
