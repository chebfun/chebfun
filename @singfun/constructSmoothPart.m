function s = constructSmoothPart(op, vscale, hscale, pref)
%CONSTRUCTSMOOTHPART   Construct the smooth part of a SINGFUN.
%   S = CONSTRUCTSMOOTHPART(OP, VSCALE, HSCALE, PREF) creates the SMOOTHFUN 
%   object S using the function handle OP, which is assumed globally smooth. 
%   VSCALE, HSCALE and preferences to the SMOOTHFUN class are passed through
%   the corresponding arguments.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

pref = chebtech.pref(pref, 'sampletest', 0);

% Call the CHEBTECH constructor:
s = smoothfun.constructor(op, vscale, hscale, pref);

end