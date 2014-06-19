function out = issing(f)
%ISSING   Test if a DELTAFUN object is built upon SINGFUN.
%   out = ISSING(F) returns logical true if the FUNPART of F is made of a 
%   SINGFUN and false otherwise.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~issing(f.funPart) )
    out = false;
else
    out = true;
end

end
