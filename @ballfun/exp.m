function g = exp(f)
%EXP   Exponential of a BALLFUN.
%   EXP(F) computes the exponential of the BALLFUN F.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @exp ); 

end
