function d = domain(f)
%DOMAIN   Returns the domain of a DELTAFUN.
%   D = DOMAIN(F) returns the domain of the FUNPAT of the DELTAFUN F. 
%
% See also MERGECOLUMNS, MERGEDELTAS, CLEANROWS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

d = get(f.funPart, 'domain');

end
