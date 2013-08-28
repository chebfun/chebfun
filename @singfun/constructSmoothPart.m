function s = constructSmoothPart(op, pref)
%CONSTRUCTSMOOTHPART Construct the smooth part of a SINGFUN.
%
%   S = CONSTRUCTSMOOTHPART(OP, PREF) creates the SMOOTHFUN object S using 
%   the function handle OP, which is assumed globally smooth. Preferences 
%   to the SMOOTHFUN class are passed through PREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.


% Construct the smooth part of the SINGFUN object.
% [TODO]: This will be replaced by a call to the SMOOTHFUN constructor

% get CHEBTECH perferences if provided
if (isfield( pref, 'chebtech') )
    pref = pref.chebtech;
else
    % get preferences from the CHEBTECH class
    pref = chebtech.pref;
end
%smoothPrefs = chebtech.pref('extrapolate', false);
vscale = [];
hscale = [];
s = chebtech.constructor(op, vscale, hscale, pref);