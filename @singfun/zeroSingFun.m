function s = zeroSingFun()
%ZEROSUNGFUN constructs the zero singfun
s = singfun;
prefs = chebtech.pref('tech', 'cheb1', 'extrapolate', false);
vscale = [];
hscale = [];
s.smoothPart = chebtech.constructor(0, vscale, hscale, prefs);
s.exponents = [0, 0];
s.isSingEnd = [0, 0];
s.singType = {'none', 'none'};

end