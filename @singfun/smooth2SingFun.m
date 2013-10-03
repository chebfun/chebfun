function f = smooth2SingFun(f)
%SMOOTH2SINGFUN   Convert a SMOOTHFUN to a SINGFUN.
%   smooth2SINGFUN(F) Converts the SMOOTHFUN object F to a SINGFUN F.
%
% See also MAKE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Take no action if f is already a SINGFUN:
if ( isa(f, 'singfun') )
    return
end

% Convert a SMOOTHFUN to a SINGFUN object:
if ( isa(f, 'smoothfun') )
    g = f;
    f = singfun.zeroSingFun();
    f.smoothPart = g; 
else
    error('CHEBFUN:SINGFUN:CONVERT2SINGFUN', 'Unknown type of function passed.');
end

end