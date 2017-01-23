function tr = trace(f)
%TRACE   Integral of a DISKFUN along its diagonal.
%   TRACE(F) is the integral of the radial slice of the DISKFUN F occuring 
%   at the angle theta = pi/4.
%
% See also DISKFUN/DIAG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    tr = [];
    return
end 

tr = sum(feval(f, pi/4, ':'));

end