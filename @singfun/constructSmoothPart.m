function s = constructSmoothPart(op, pref)
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