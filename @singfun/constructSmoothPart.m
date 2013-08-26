function s = constructSmoothPart(op, pref)
% Construct the smooth part of the SINGFUN object.
% [TODO]: This will be replaced by a call to the SMOOTHFUN constructor

pref = pref.chebtech;
%smoothPrefs = chebtech.pref('extrapolate', false);
vscale = [];
hscale = [];
s = chebtech.constructor(op, vscale, hscale, pref);