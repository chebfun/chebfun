function f = fracInt(f, mu)
%FRACINT  Fractional integral of an UNBNDFUN. 
%   Fractional integrals on unbounded domain are not supported.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:UNBNDFUN:fracInt:noSupport', ...
    'Fractional integrals on unbounded domains are not yet supported.')

end
