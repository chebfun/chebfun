function s = zeroDeltaFun(domain)
%ZEROSINGFUN   Constructs the zero DELTAFUN on DOMAIN. The output DELTAFUN 
%   object has a trivial smooth part with no delta functions.
%
% See also DELTAFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Create a zero fun on the domain:
f = fun.constructor(0, domain);

% Create a zero DELTAFUN object:
s = deltafun();
s.funPart = f;
s.impulses = [];
s.location = [];
end
