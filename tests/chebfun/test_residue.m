function pass = test_residue(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

f = chebfun(@(x) (x-1.1).*(x.^2+1).*(x-10i), pref);
g = chebfun(@(x) x.^5, pref);

[r, p, k] = residue(g, f);
[G, F] = residue(r, p, k);

tol = 100*eps;

pexact = [10i ; 1.1 ; -1i ; 1i];

pass(1) = norm(sort(real(p)) - sort(real(pexact))) + ...
          norm(sort(imag(p)) - sort(imag(pexact))) < 10*tol;

pass(2) = normest(g.*F - G.*f)./vscale(g.*F) < 10*tol*vscale(k);

pass(3) = normest(k - chebfun(@(x) x + 10i + 1.1)) < 10*tol;

% Check syntax which substitutes 0 for empty third argument.
[B, A] = residue([1 1], [1 -1], chebfun());
pass(4) = (normest(B - chebfun(@(x) 2*x)) < tol) && ...
    (normest(A - chebfun(@(x) x.^2 - 1)) < tol);

end
