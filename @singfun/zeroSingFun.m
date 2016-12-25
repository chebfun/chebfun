function s = zeroSingFun()
%ZEROSINGFUN   Constructs the zero SINGFUN. The output SINGFUN object has a
%   zero smooth part with no singularities at any end points.
%
% See also SINGFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Create an empty SINGFUN object:
s = singfun();

% Create a zero smooth part:
s.smoothPart = singfun.constructSmoothPart(0, [], []);

% No singularities at any end points:
s.exponents = [0, 0];

end
