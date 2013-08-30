function s = zeroSingFun()
%ZEROSUNGFUN   Constructs the zero SINGFUN. The output SINGFUN object has a
%   zero smooth part with no singularities at any end points.
%
% See also SINGFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Create an empty SINGFUN object:
s = singfun;

% Create a zero smooth part:
pref = chebtech.pref('tech', 'cheb1', 'extrapolate', false);
vscale = [];
hscale = [];
pref = smoothfun.pref;
s.smoothPart = smoothfun.constructor(@(x) 0*x, vscale, hscale, pref);

% No singularities at any end points:
s.exponents = [0, 0];
s.singType = {'none', 'none'};

end
