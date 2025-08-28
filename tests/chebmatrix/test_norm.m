function pass = test_norm(pref)
% HM, 30 Apr 2014

%% 
% Obtain preferences.
if ( nargin == 0 )
    pref = chebfunpref();
end

%% 
% The entries of A are only CHEBFUN or DOUBLE.
f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) cos(x), [-1 -0.5 0 0.5 1], pref);
h = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);

A = [ f ; g ; h ];

chebA = chebfun(A);
tol = 2*eps*vscale(chebA);

norm2 = 2.372100421113536830;
normInf = 2.718281828459046;

err(1) = abs(norm(A) - norm2);       
err(2) = abs(norm(A, 'fro') - norm2);   
err(3) = abs(norm(A, 2) - norm2);
err(4) = abs(norm(A, inf) - normInf);
err(5) = abs(norm(A, 'inf') - normInf);

pass = err < tol;
    
%%
% [TODO]: Add tests for norms of operators (inf x inf blocks).

end
