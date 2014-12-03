function out = isTrigTech(f)
%ISTRIGTECH   Test if a CHEBFUN object is built upon TRIGTECH.
%   out = ISTRIGTECH(F) returns logical true if F has at least one FUN
%   which is made of TRIGTECH and false otherwise.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = any(cellfun(@(f) isTrigTech(f), f.funs));

end
