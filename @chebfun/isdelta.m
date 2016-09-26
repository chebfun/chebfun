function out = isdelta(f)
%ISDELTA   Test if a CHEBFUN object is built upon DELTAFUN.
%   out = ISDELTA(F) returns logical true if F has at least one FUN which is 
%   made of DELTAFUN and false otherwise.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = any(cellfun(@(f) isdelta(f), f.funs));

end
