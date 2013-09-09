function s = plus(f, g)
%+   Addition of DELTAFUN objects.
%   F + G adds F and G, where F and G may be DELTAFUN objects or scalars.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If one of the arguments is empty:
if ( isempty(f) || isempty(g) )
    % Create an empty singfun and return:
    s = deltafun;
    return;
end

% One of the arguments i.e. f or g is necessarily a SINGFUN object. Otherwise, 
% this overloaded plus would not have been called.

% If one of the arguments is a double:
if ( isa(f, 'double') )
    aDouble = f;
    f = singfun.zeroSingFun();
    f.smoothPart = singfun.constructSmoothPart(aDouble, []);
elseif ( isa(g, 'double') )
    aDouble = g;
    g = singfun.zeroSingFun();
    g.smoothPart = singfun.constructSmoothPart(aDouble, []);
end


end