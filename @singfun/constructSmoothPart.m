function s = constructSmoothPart(op, pref)
%CONSTRUCTSMOOTHPART   Construct the smooth part of a SINGFUN.
%   S = CONSTRUCTSMOOTHPART(OP, PREF) creates the SMOOTHFUN object S using the
%   function handle OP, which is assumed globally smooth. Preferences to the
%   SMOOTHFUN class are passed through PREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: This will be replaced by a call to the SMOOTHFUN constructor

if ( isfield( pref, 'chebtech') )
    % Get CHEBTECH preferences if provided:
    pref = pref.chebtech;
else
    % Get default preferences from the CHEBTECH class:
    pref = chebtech.pref();
end

% Call the chebtech constructor:
s = chebtech.constructor(op, [], [], pref);

end
