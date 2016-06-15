function L = laplacian(F)
%LAPLACIAN   Vector Laplacian of a CHEBFUN3V object.
%   LAPLACIAN(F) returns a CHEBFUN3V representing the vector Laplacian of 
%   F. The output L is the vector field of the scalar Laplacian applied to
%   the individual components.
%
% See also CHEBFUN3V/LAP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) )
    L = chebfun3v;
    return
end

diff1 = diff(F, 2, 1);
diff2 = diff(F, 2, 2); 
diff3 = diff(F, 2, 3);

dom = diff1.components{1}.domain;

% laplacian = F_xx + F_yy + F_zz
L = diff1; % Initialize the output

% Compute 1st component of the output:
vscales1 = [vscale(diff1.components{1}) + vscale(diff2.components{1}), ...
    vscale(diff3.components{1})];
% Developer Note: Instead of calling e.g. L1 = diff1 + diff2 + diff3; which
% needs to call the constructor twice, we use the following to call it just
% once. See CHEBFUN3/PLUS for more details.
m = 51; % size of sampling grid
L1Vals = sample(diff1.components{1}, m, m, m) + sample(...
    diff2.components{1}, m, m, m) + sample(diff3.components{1}, m, m, m);
L1Vscale = max(abs(L1Vals(:)));
kappa1 = sum(vscales1)/L1Vscale;
pref = chebfunpref(); prefStruct = pref.cheb3Prefs;
eps = prefStruct.chebfun3eps;
tol = eps*kappa1;
L.components{1} = chebfun3(@(x,y,z) feval(diff1.components{1}, x, y, z) + ...
    feval(diff2.components{1}, x, y, z) + feval(diff3.components{1}, ...
    x, y, z) , dom, 'eps', tol);

% Compute 2nd component of the output:
vscales2 = [vscale(diff1.components{2}) + vscale(diff2.components{2}), ...
    vscale(diff3.components{2})];
L2Vals = sample(diff1.components{2}, m, m, m) + sample(...
    diff2.components{2}, m, m, m) + sample(diff3.components{2}, m, m, m);
L2Vscale = max(abs(L2Vals(:)));
kappa2 = sum(vscales2)/L2Vscale;
tol = eps*kappa2;
L.components{2} = chebfun3(@(x,y,z) feval(diff1.components{2}, x, y, z) + ...
    feval(diff2.components{2}, x, y, z) + feval(diff3.components{2}, ...
    x, y, z) , dom, 'eps', tol);

% Compute 3rd component of the output:
vscales3 = [vscale(diff1.components{3}) + vscale(diff2.components{3}), ...
    vscale(diff3.components{3})];
L3Vals = sample(diff1.components{3}, m, m, m) + sample(...
    diff2.components{3}, m, m, m) + sample(diff3.components{3}, m, m, m);
L3Vscale = max(abs(L3Vals(:)));
kappa3 = sum(vscales3)/L3Vscale;
tol = eps*kappa3;
L.components{3} = chebfun3(@(x,y,z) feval(diff1.components{3}, x, y, z) + ...
    feval(diff2.components{3}, x, y, z) + feval(diff3.components{3}, ...
    x, y, z) , dom, 'eps', tol);

end