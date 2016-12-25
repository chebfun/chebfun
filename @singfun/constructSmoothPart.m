function s = constructSmoothPart(op, data, pref)
%CONSTRUCTSMOOTHPART   Construct the smooth part of a SINGFUN.
%   S = CONSTRUCTSMOOTHPART(OP, DATA, PREF) creates the SMOOTHFUN object S
%   using the function handle OP, which is assumed to be globally smooth.  Data
%   and preferences for the SMOOTHFUN class can be passed through the
%   corresponding arguments.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(pref) )
    pref = chebfunpref();
end

% Call the SMOOTHFUN constructor:
s = smoothfun.constructor(op, data, pref);

end
