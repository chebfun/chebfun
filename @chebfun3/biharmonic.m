function B = biharmonic(f)
%BIHARMONIC   Biharmonic operator applied to a CHEBFUN3.
%   B = BIHARMONIC(F) returns a CHEBFUN3 object B representing the 
%   biharmonic operator applied to F.
%
% See also CHEBFUN3/BIHARM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% biharmonic(f) = f_xxxx + f_yyyy + f_zzzz + 2*f_xxyy + 2*f_xxzz + 2*f_yyzz.

diff4x = diff(f, 4, 1); 
diff4y = diff(f, 4, 2); 
diff4z = diff(f, 4, 3);
diff2x2y = 2*diff(diff(f, 2, 1), 2, 2);
diff2x2z = 2*diff(diff(f, 2, 1), 2, 3);
diff2y2z = 2*diff(diff(f, 2, 2), 2, 3);

% See CHEBFUN3/PLUS for more details.
vscales = [vscale(diff4x), vscale(diff4y), vscale(diff4z), ...
    vscale(diff2x2y), vscale(diff2x2z), vscale(diff2y2z)];
m = 51; % size of sampling grid
BVals = sample(diff4x, m, m, m) + sample(diff4y, m, m, m) + ...
    sample(diff4z, m, m, m) + sample(diff2x2y, m, m, m) + ...
    sample(diff2x2z, m, m, m) + sample(diff2y2z, m, m, m);

BVscale = max(abs(BVals(:)));
kappa = sum(vscales)/BVscale;
pref = chebfunpref().cheb3Prefs;
eps = pref.chebfun3eps;
tol = eps*kappa;
B = chebfun3(@(x,y,z) feval(diff4x, x, y, z) + feval(diff4y, x, y, z) + ...
    feval(diff4z, x, y, z) + feval(diff2x2y, x, y, z) + ...
    feval(diff2x2z, x, y, z) + feval(diff2y2z, x, y, z), f.domain, 'eps', tol);


end