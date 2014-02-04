function d = domain(f)
%DOMAIN   Returns the domain of a DELTAFUN.
%   D = DOMAIN(F) returns the domain of the FUNPAT of the DELTAFUN F. 
%
% See also MERGECOLUMNS, MERGEIMPULSES, CLEANROWS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

d = f.funPart.domain;

end