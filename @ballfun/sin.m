function g = sin(f)
%SIN  Sine of a BALLFUN.
%   SIN(F) computes the sine of the BALLFUN F.
%
% See also SINH, COS.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @sin ); 

end
