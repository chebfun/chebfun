function out = issing(f)
%ISSING   Test if a CHEBFUN object is built upon SINGFUN.
%   out = ISSING(F) returns logical true if F has at least one FUN which is 
%   made of SINGFUN and false otherwise.

% [TODO]: Should the name of this function go camel case? So far, all is*
% function are named without camel case.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = any(cellfun(@(f) issing(f), f.funs));

end
