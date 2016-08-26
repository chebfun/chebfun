function pass = test_eig( pref ) 
% Test for EIG command of Chebfun2.

if ( nargin == 0 ) 
    pref = chebfunpref;
end 
tol = 1e4*pref.cheb2Prefs.chebfun2eps;

% Decomposition on a square domain.
f = cheb.gallery2('challenge');
[V, D] = eig(f);

% Is it correct?
pass(1) = norm(f * V - V * D) < tol;

%%
% The following should give an error message as the domain is not square.
f = chebfun2(@(x,y,z) x+y, [-1 1 -2 2]);
try
    d = eig(f);
    pass(2) = false;
catch ME 
    pass(2) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN2:eig:domainerr');
end

end