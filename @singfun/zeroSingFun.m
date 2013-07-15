function s = zeroSingFun()
%ZEROSUNGFUN   Constructs the zero SINGFUN. The output SINGFUN object has a
%   zero smooth part with no singularities at any end points.
%
% See also SINGFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%
% create an empty SINGFUN object.
s = singfun;

%%
% creat a zero smooth part
prefs = chebtech.pref('tech', 'cheb1', 'extrapolate', false);
vscale = [];
hscale = [];
s.smoothPart = chebtech.constructor(0, vscale, hscale, prefs);

%%
% No singularities at any end points.
s.exponents = [0, 0];
s.isSingEnd = [0, 0];
s.singType = {'none', 'none'};

end