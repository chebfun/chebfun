function out = isunbnd(f)
%ISUNBND   Test if a FUN object is defined on an unbounded domain.
%   out = ISUNBND(F) returns logical true if F is defined on an unbounded
%   domain, singly unbounded or doubly unbounded.

% [TODO]: Should the name of this function go camel case? So far, all is*
% function are named without camel case.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = any( isinf(f.domain) );

end
