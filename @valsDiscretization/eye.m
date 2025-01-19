function I = eye(disc)
%EYE   Identity operator for VALSDISCRETIZATION discretization.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = disc.dimension;
pref = cheboppref;
if ( pref.sparse )
    I = speye(sum(n));
else
    I = eye(sum(n));
end

end
