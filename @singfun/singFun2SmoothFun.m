function f = singFun2SmoothFun(f)
%SINGFUN2SMOOTHFUN   Convert a SINGFUN to a SMOOTHFUN, if possible.
%   SINGFUN2SMOOTHFUN(F) converts the SINGFUN object F to a SMOOTHFUN F, if the 
%   exponents are negligible.
%
% See also MAKE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( issmooth(f) )
    f = f.smoothPart; 
end

end
