function out = issing(f)
%ISSING   Test if a FUN object is built upon SINGFUN.
%   out = ISSING(F) returns logical true if its ONEFUN is a SINGFUN and false 
%   otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

out = isa(f.onefun, 'singfun');

end