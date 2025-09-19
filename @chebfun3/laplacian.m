function L = laplacian(f)
%LAPLACIAN   Laplacian of a CHEBFUN3.
%   L = LAPLACIAN(F) returns a CHEBFUN3 representing the Laplacian of F.
%
% See also CHEBFUN3/LAP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

diff1 = diff(f, 2, 1); 
diff2 = diff(f, 2, 2); 
diff3 = diff(f, 2, 3);

vscales = vscale(diff1) + vscale(diff2) + vscale(diff3);
% Developer Note: Instead of calling 
% L = diff1 + diff2 + diff3; which needs to call the constructor twice, we
% use the following to call it just once. See CHEBFUN3/PLUS for more
% details.

if vscales == 0
    L = chebfun3(0);
else
    % Find out a proper tolerance for the plus operation
    m = 51; % size of sampling grid
    LVals = sample(diff1, m, m, m) + sample(diff2, m, m, m) + ...
        sample(diff3, m, m, m);
    LVscale = max(abs(LVals(:)));
    kappa = vscales/LVscale;
    pref = chebfunpref();
    eps = pref.cheb3Prefs.chebfun3eps;
    tol = eps*kappa;

    % Compute the laplacian:
    L = chebfun3(@(x,y,z) feval(diff1, x, y, z) + feval(diff2, x, y, z) + ...
        feval(diff3, x, y, z) , f.domain, 'eps', tol);
end

end