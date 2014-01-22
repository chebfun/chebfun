function f = singFun2SmoothFun(f)
%SINGFUN2SMOOTHFUN   Convert a SINGFUN to a SMOOTHFUN.
%   SINGFUN2SMOOTHFUN(F) Converts the SINGFUN object F to a SMOOTHFUN F.
%
% See also MAKE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Take no action if f is already a SMOOTHFUN:
if ( isa(f, 'smoothfun') )
    return
end

% Convert a SINGFUN to a SMOOTHFUN object if EXPONENTS of F are trivial:
if ( isa(f, 'singfun') && ~any(f.exponents) )
    f = f.smoothPart; 
end

end