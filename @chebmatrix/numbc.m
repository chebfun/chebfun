function k = numbc(L)
% NUMBC(L) returns the number of constraints attached to L.
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
k = length(L.constraints);
end
