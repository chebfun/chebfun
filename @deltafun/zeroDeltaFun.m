function s = zeroDeltaFun(domain)
%ZEROSINGFUN   Constructs the zero DELTAFUN on DOMAIN. The output DELTAFUN 
%   object has a trivial smooth part with no delta functions.
%
% See also DELTAFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Create a zero smooth part:
f = chebfun(0, domain);

% Create zero DELTAFUN object:
% NOTE: To avoid constructor nuances, we are placing a 0-mag deltafunction at
% the mid-point of the domain.
a = domain(1);
b = domain(2);
s = deltafun(0, (a+b)/2, f, domain);
end
