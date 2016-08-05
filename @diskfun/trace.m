function tr = trace(f)
% TRACE integral of a DISKFUN along its diagonal 
%
% TRACE(f) is the integral of the radial slice of the disk
% occuring at the angle theta=pi/4. 
% 
% See also DIAG.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) ) 
    tr = [];
    return
end 

tr = sum(feval(f, pi/4, ':'));

end
