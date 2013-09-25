function s = constructSmoothPart(op, vscale, hscale, pref)
%CONSTRUCTSMOOTHPART   Construct the smooth part of a SINGFUN.
%   S = CONSTRUCTSMOOTHPART(OP, PREF) creates the SMOOTHFUN object S using the
%   function handle OP, which is assumed globally smooth. Preferences of the
%   SMOOTHFUN class and underlying CHEBTECH preferences are passed through PREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% If none of SMOOTHFUN and CHEBTECH pref is supplied:
if ( ( nargin < 4 ) || (~isfield(pref, 'smoothfun') && ~isfield(pref, 'chebtech')) )
    pref = smoothfun.pref;
    pref = chebtech.pref(pref);
end

% If only CHEBTECH pref is not supplied:
if ( isfield(pref, 'smoothfun') && ~isfield(pref, 'chebtech') )
    % Get default CHEBTECH preferences if not provided:
    pref = chebtech.pref(pref);
end

% If only SMOOTHFUN pref is not supplied:
if ( isfield(pref, 'chebtech') && ~isfield(pref, 'smoothfun') )
    % Get default SMOOTHFUN preferences if not provided:
    pref = smoothfun.pref(pref);
end 

% For time-being, SAMPLETEST is turned off:
pref = chebtech.pref(pref, 'sampletest', 0);

% Though both 'CHEBTECH1' and 'CHEBTECH2' work well for singfun, in long term we
% are aiming at using first kind Chebyshev points:
pref = chebtech.pref(pref, 'tech', 'chebtech1');

% Call the CHEBTECH constructor:
s = smoothfun.constructor(op, vscale, hscale, pref);

end