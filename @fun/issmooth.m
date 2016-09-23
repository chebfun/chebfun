function out = issmooth(f)
%ISSMOOTH   Test if a FUN object is built upon SMOOTHFUN.
%   OUT = ISSMOOTH(F) returns logical true if its ONEFUN is a SMOOTHFUN and 
%   false otherwise.

% [TODO]: Should the name of this function go camel case? So far, all is*
% function are named without camel case.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = isa(f.onefun, 'smoothfun');

end
