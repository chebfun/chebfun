function pass = test_domainvolume(pref)
% Test chebfun3/domainvolume

if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

% Default domain: 
ff = @(x,y,z) cos(x+y+z);
f = chebfun3(ff);
pass(1) = domainvolume(f) == 8;

% A different domain: 
dom = [-1 2 -2 1 -3 0];
f = chebfun3(ff, dom);
pass(2) = domainvolume(f) == 27;

% Another domain: 
dom = [-pi pi -pi pi -pi pi];
f = chebfun3(ff, 'trig', dom);
pass(3) = abs(domainvolume(f) - (2*pi)^3) < tol;

end 