function out = angles(n)
%ANGLES   Return the angles of the Chebyshev points of 2nd kind in [-1, 1].
%   CHEBTECH2.ANGLES(N) returns ACOS(X), where X are the N Chebyshev points of
%   the 2nd kind in [-1, 1].
%
% See also POINTS, CHEBPTS, LENGTH.
%
% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

if ( n == 0 )
    out = [];
    return
elseif ( n == 1 )
    out = pi/2;
    return 
end

m = n - 1;
out = (m:-1:0).'*pi/m;

end
