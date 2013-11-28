% Test file for chebfun/innerProduct.m

function pass = test_innerProduct(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Set a domain
dom = [-2 7];

%% Integration with singfun

pow1 = -0.3;
pow2 = -0.5;
op1 = @(x) (x - dom(2)).^pow1.*sin(x);
op2 = @(x) (x - dom(2)).^pow2.*cos(3*x);
pref.singPrefs.exponents = [0 pow1];
f = chebfun(op1, dom, pref);
pref.singPrefs.exponents = [0 pow2];
g = chebfun(op2, dom, pref);
I = innerProduct(f,g);
I_exact = -0.65182492763883119+0.47357853074362785i;
pass(1) = ( abs(I-I_exact) < 1e1*max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

end